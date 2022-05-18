import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import single_reaction_deletion
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis
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
                    ""
                    self.dict_biocyc2metanetx[line[0]] = metanetx

    def add_rxn(self, list_of_rxns):
        list_of_rxns = list_of_rxns.split("\n")
        self.model_copy = self.model.copy()
        self.dict_metabolite2 = copy.deepcopy(self.dict_metabolite2biocyc)
        isLastRxn = False

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

    def add_metabolite(self, biocyc_unique_id):
        metabolite = Metabolite(biocyc_unique_id)
        self.dict_metabolite2[biocyc_unique_id] = metabolite
        self.model_copy.add_metabolites(metabolite)
        return metabolite

    def run(self, list_of_pathways):
        for pathway in list_of_pathways:
            self.add_rxn(pathway)
            linear_reaction_coefficients(self.model)
            theoretical_yield = self.model_copy.optimize().objective_value

            nadh_summary = self.model_copy.metabolites.nadh_c.summary()
            nadh_consumption = nadh_summary.consuming_flux["flux"].sum()

            nadph_summary = self.model_copy.metabolites.nadph_c.summary()
            nadph_consumption = nadph_summary.consuming_flux["flux"].sum()

            atp_summary = self.model_copy.metabolites.atp_c.summary()
            atp_consumption = atp_summary.consuming_flux["flux"].sum()

            print(str(theoretical_yield) + "\t" + str(atp_consumption)+ "\t" + str(nadh_consumption) + "\t" + str(nadph_consumption))
            self.model_copy = None


def runner():
    model_path = sys.argv[1]
    list_of_paths = sys.argv[2]
    list_of_paths = list_of_paths.split("//")
    converter = CobraConverter(model_path)
    converter.run(list_of_paths)


runner()

if __name__ == "__main__":
    runner(
        "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/iJO1366.xml",
        "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n"
        "RXN-14985\tCPD-3618 PROTON  --> CARBON-DIOXIDE BUTANAL\n"
        "RXN-14986\tNAD CPD-1130  --> NADH CPD-3618 CARBON-DIOXIDE\n"
        "3-ETHYLMALATE-SYNTHASE-RXN\tWATER GLYOX BUTYRYL-COA  --> CO-A CPD-1130 PROTON\n"
        "RXN-11726\tNADPH CROTONYL-COA PROTON  --> NADP BUTYRYL-COA\n"
        "3-HYDROXBUTYRYL-COA-DEHYDRATASE-RXN\tCPD-650  --> WATER CROTONYL-COA\n"
        "RXN-5901\tACETOACETYL-COA NADPH PROTON  --> NADP CPD-650\n" + "//" +
        "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n" 
        "RXN-14985\tCPD-3618 PROTON  --> CARBON-DIOXIDE BUTANAL\n"
        "RXN-14986\tNAD CPD-1130  --> NADH CPD-3618 CARBON-DIOXIDE\n"
        "3-ETHYLMALATE-SYNTHASE-RXN\tWATER GLYOX BUTYRYL-COA  --> CO-A CPD-1130 PROTON\n"
        "2.8.3.9-RXN\tACETOACETYL-COA BUTYRIC_ACID  --> 3-KETOBUTYRATE BUTYRYL-COA\n"
        "RXN-11662\tNAD S-3-HYDROXYBUTANOYL-COA  --> NADH ACETOACETYL-COA PROTON\n"
    )
