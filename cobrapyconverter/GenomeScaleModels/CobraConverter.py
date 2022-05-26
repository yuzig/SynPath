import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import single_reaction_deletion
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion
import copy
import sys

class CobraConverter:
    dict_metabolite2biocyc = {}
    dict_metabolite2metanetx = {}
    dict_biocyc2metanetx = {}
    dict_rxn = {}
    dict_metabolite2 = {}
    model = None
    model_copy = None
    last_reactions = []
    growth_rate = 0
    atpm = 0

    def __init__(self, model_file_path):

        self.model = cobra.io.read_sbml_model(model_file_path)
        for r in self.model.reactions:
            if r.annotation.get('biocyc') is not None and type(r.annotation.get('biocyc')) == str:
                self.dict_rxn[r.annotation.get('biocyc')] = r
            if r.annotation.get('biocyc') is not None and type(r.annotation.get('biocyc')) == list:
                for biocyc_id in r.annotation.get('biocyc'):
                    self.dict_rxn[biocyc_id] = r

        for c in self.model.metabolites:
            if c.annotation.get('biocyc') is not None and type(c.annotation.get('biocyc')) == str:
                self.dict_metabolite2biocyc[c.annotation.get('biocyc')] = c
            if c.annotation.get('biocyc') is not None and type(c.annotation.get('biocyc')) == list:
                for biocyc_id in c.annotation.get('biocyc'):
                    self.dict_metabolite2biocyc[biocyc_id] = c
            if c.annotation.get('metanetx.chemical') is not None:
                self.dict_metabolite2metanetx[c.annotation.get('metanetx.chemical')] = c
        self.dict_metabolite2biocyc["META:NADH-P-OR-NOP"] = self.model.metabolites.get_by_id('nadh_c')
        self.dict_metabolite2biocyc["META:NAD-P-OR-NOP"] = self.model.metabolites.get_by_id('nad_c')
        self.dict_metabolite2biocyc["META:Acceptor"] = self.model.metabolites.get_by_id('fad_c')
        self.dict_metabolite2biocyc["META:Donor-H2"] = self.model.metabolites.get_by_id('fadh2_c')

        with open(
                '/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/biocyc2metanetx.txt') as f:
            lines = f.readlines()
            for line in lines:
                line = line.replace("\n", "")
                line = line.split("\t")
                if line[1] == "null":
                    continue
                else:
                    metanetx = line[1].replace('"', '')
                    self.dict_biocyc2metanetx[line[0]] = metanetx
        with open(
                '/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/biocyc2bigg.txt') as f:
            lines = f.readlines()
            for line in lines:
                line = line.replace("\n", "")
                line = line.split("\t")
                if line[1] == "null":
                    continue
                else:
                    bigg_id = line[1].replace('"', '')
                    bigg_id = bigg_id + "_c"
                    try:
                        self.dict_metabolite2biocyc[line[0]] = self.model.metabolites.get_by_id(bigg_id)
                    except:
                        continue

        linear_reaction_coefficients(self.model)
        solution = self.model.optimize()
        self.growth_rate = solution.objective_value

    def add_rxn(self, list_of_rxns):
        list_of_rxns = list_of_rxns.split("\n")
        self.model_copy = self.model.copy()
        self.model_copy.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M").lower_bound = 0.8 * self.growth_rate

        self.dict_metabolite2 = copy.deepcopy(self.dict_metabolite2biocyc)
        isLastRxn = False
        for i in range(len(list_of_rxns)):
            if list_of_rxns[i] == "":
                list_of_rxns.pop(i)
        possible = False

        for i in range(1, len(list_of_rxns)):
            if self.is_native(list_of_rxns[len(list_of_rxns) - i]):
                possible = True
                break

        if not possible:
            return False

        for i in range(len(list_of_rxns)):
            if i == 0:
                isLastRxn = True
            r = list_of_rxns[i]

            if r == '' or r == ' ':
                continue
            tabs = r.split("\t")
            t = tabs[0].replace(" ", "")
            reaction = Reaction(t)

            if "META:" + tabs[0] in self.dict_rxn or tabs[0] == '':
                if isLastRxn:
                    self.model_copy.objective = self.dict_rxn.get("META:" + tabs[0]).id
                    self.last_reactions.append(reaction)
                    isLastRxn = False
                continue

            self.model_copy.add_reactions({reaction})

            if isLastRxn:
                self.model_copy.objective = reaction.id
                self.last_reactions.append(reaction)

            metabolites = tabs[1].split(" --> ")
            reactants = metabolites[0].split(' ')
            products = metabolites[1].split(' ')
            for c in reactants:
                if c == '':
                    continue
                if c == 'Red-NADPH-Hemoprotein-Reductases' or c == 'Ox-NADPH-Hemoprotein-Reductases':
                    continue
                if "META:" + c in self.dict_metabolite2:
                    c = self.dict_metabolite2.get("META:" + c)
                    reaction.add_metabolites({c.id: -1})
                    continue

                if self.dict_biocyc2metanetx.get(c) in self.dict_metabolite2metanetx:
                    c = self.dict_metabolite2metanetx.get(self.dict_biocyc2metanetx.get(c))
                    reaction.add_metabolites({c.id: -1})
                else:
                    metabolite = self.add_metabolite("META:" + c)
                    reaction.add_metabolites({metabolite.id: -1})

            for c in products:
                if c == '':
                    continue
                if "META:" + c in self.dict_metabolite2:
                    c = self.dict_metabolite2.get("META:" + c)
                    reaction.add_metabolites({c.id: 1})
                    continue
                if self.dict_biocyc2metanetx.get(c) in self.dict_metabolite2metanetx:
                    c = self.dict_metabolite2metanetx.get(self.dict_biocyc2metanetx.get(c))
                    reaction.add_metabolites({c.id: 1})
                else:
                    metabolite = self.add_metabolite("META:" + c)
                    reaction.add_metabolites({metabolite.id: 1})
                    if isLastRxn:
                        self.model_copy.add_boundary(self.model_copy.metabolites.get_by_id("META:" + c), type="demand")

            isLastRxn = False
        return True

    def is_native(self, reaction_str):
        metabolites = reaction_str.split(" --> ")
        metabolites = metabolites[0].split("\t")
        reactants = metabolites[1].split(' ')
        for c in reactants:
            if c == '':
                continue
            if not self.dict_biocyc2metanetx.get(c) in self.dict_metabolite2metanetx and not "META:" + c in self.dict_metabolite2:
                return False
        return True

    def add_metabolite(self, biocyc_unique_id):
        metabolite = Metabolite(biocyc_unique_id)
        self.dict_metabolite2[biocyc_unique_id] = metabolite
        self.model_copy.add_metabolites(metabolite)
        return metabolite

    def run(self, list_of_pathways):
        idx = 0
        for pathway in list_of_pathways:
            if not self.add_rxn(pathway):
                idx = idx + 1
                continue

            # Theoretical yield calculation
            linear_reaction_coefficients(self.model)
            theoretical_yield = self.model_copy.optimize().objective_value
            if theoretical_yield == 0.0:
                idx = idx + 1
                continue

            # deletion_results = single_reaction_deletion(self.model_copy)
            #
            # i = 0
            # deletion_set = set()
            # while i <= 3:
            #     idx = deletion_results["growth"].idxmax()
            #     max_gene = deletion_results.loc[idx, "ids"]
            #     max_gene1 = "".join(max_gene)
            #     max_gene1 = self.model_copy.reactions.get_by_id(max_gene1)
            #     deletion_set.add(max_gene1)
            #     deletion_results.drop(idx,axis=0,inplace=True)
            #     i = i + 1
            #
            # self.model_copy.remove_reactions(deletion_set)
            # eng_yield = self.model_copy.optimize().objective_value
            #
            # # cofactor calculation
            # nadh_summary = self.model_copy.metabolites.nadh_c.summary()
            # nadh_consumption = nadh_summary.consuming_flux["flux"].sum()
            #
            # # fmnh2_summary = self.model_copy.metabolites.fmnh2_c.summary()
            # # fmnh2_consumption = fmnh2_summary.consuming_flux["flux"].sum()
            #
            # # fadh2_summary = self.model_copy.metabolites.fadh2_c.summary()
            # # fadh2_consumption = fadh2_summary.consuming_flux["flux"].sum()
            #
            # nadph_summary = self.model_copy.metabolites.nadph_c.summary()
            # nadph_consumption = nadph_summary.consuming_flux["flux"].sum()
            #
            # atp_summary = self.model_copy.metabolites.atp_c.summary()
            # atp_consumption = atp_summary.consuming_flux["flux"].sum()
            # #
            # # FVA span calculation
            # # self.model_copy.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M").lower_bound = 0.5 * self.growth_rate
            # # FVA = flux_variability_analysis(self.model_copy)
            # # dif = FVA["maximum"] - FVA["minimum"]
            # # fva = dif.sum()

            atp_consumption = 0
            nadh_consumption = 0
            nadph_consumption = 0
            eng_yield = 0

            print(str(theoretical_yield)
                  + "\t" + str(atp_consumption)
                  + "\t" + str(nadh_consumption)
                  + "\t" + str(nadph_consumption)
                  + "\t" + str(eng_yield)
                  + "\t" + str(idx))

            self.model_copy = None
            idx = idx + 1


def runner():
    model_path = sys.argv[1]
    list_of_paths = sys.argv[2]
    list_of_paths = list_of_paths.split("//")
    converter = CobraConverter(model_path)
    converter.run(list_of_paths)

#
runner()

if __name__ == "__main__":
    runner(
        "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/iML1515.xml",
        "RXN-13432\tRed-NADPH-Hemoprotein-Reductases CPD-7554 OXYGEN-MOLECULE  --> CPD-13248 PROTON Ox-NADPH-Hemoprotein-Reductases WATER\n"
        "RXN-8046\tFARNESYL-PP  --> PPI CPD-7554\n"
        "//RXN-12314\tRed-NADPH-Hemoprotein-Reductases CPD-7557 OXYGEN-MOLECULE  --> Ox-NADPH-Hemoprotein-Reductases CPD-13248 PROTON WATER\n"
        "RXN-8051\tRed-NADPH-Hemoprotein-Reductases CPD-7556 OXYGEN-MOLECULE  --> Ox-NADPH-Hemoprotein-Reductases CPD-7557 WATER\n"
        "RXN-8050\tRed-NADPH-Hemoprotein-Reductases CPD-7554 OXYGEN-MOLECULE  --> CPD-7556 Ox-NADPH-Hemoprotein-Reductases WATER\n"
        "RXN-8046\tFARNESYL-PP  --> PPI CPD-7554\n"
    )