import sys

import cobra
from IPython.core.display import display
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis, single_reaction_deletion, moma
import copy
import numpy as np
from numpy import loadtxt


class CobraConverter:
    dict_metabolite2biocyc = {}
    dict_metabolite2metanetx = {}
    dict_biocyc2metanetx = {}
    dict_rxn = {}
    dict_metabolite2 = {}
    dict_rxn_coefficients = {}
    model = model_copy = None
    growth_rate = 0
    ATPS4r = THD2 = NADH16 = biomass_rxn = None
    cofactor_recyc_val = []

    def __init__(self, model_file_path):

        self.model = cobra.io.read_sbml_model(model_file_path)
        for r in self.model.reactions:
            if r.annotation.get('biocyc') is not None and type(r.annotation.get('biocyc')) == str:
                self.dict_rxn[r.annotation.get('biocyc')] = r
            if r.annotation.get('biocyc') is not None and type(r.annotation.get('biocyc')) == list:
                for biocyc_id in r.annotation.get('biocyc'):
                    self.dict_rxn[biocyc_id] = r
            if "BIOMASS" in r.id and not "WT" in r.id:
                self.biomass_rxn = r.id
            if 'ATPS4r' in r.id:
                self.ATPS4r = r.id
            if 'THD2' in r.id:
                self.THD2 = r.id
            if 'NADH16' in r.id:
                self.NADH16 = r.id

        for c in self.model.metabolites:
            if c.annotation.get('biocyc') is not None and type(c.annotation.get('biocyc')) == str:
                self.dict_metabolite2biocyc[c.annotation.get('biocyc')] = c
            if c.annotation.get('biocyc') is not None and type(c.annotation.get('biocyc')) == list:
                for biocyc_id in c.annotation.get('biocyc'):
                    self.dict_metabolite2biocyc[biocyc_id] = c
            if c.annotation.get('metanetx.chemical') is not None:
                self.dict_metabolite2metanetx[c.annotation.get('metanetx.chemical')] = c
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

        with open('/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/reactions.txt') as f:
            lines = f.readlines()
            coefficient_dict = {}
            reaction = None
            chem_list = []
            coefficient_list = []
            for line in lines:
                if line.startswith('UNIQUE-ID'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    reaction = line[1]
                    continue
                if line.startswith('^COEFFICIENT'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    try:
                        coefficient_list[-1] = int(line[1]) * coefficient_list[-1]
                    except:
                        continue
                    continue
                if line.startswith('LEFT'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    chem_list.append(line[1])
                    coefficient_list.append(1)
                    continue
                if line.startswith('RIGHT'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    chem_list.append(line[1])
                    coefficient_list.append(1)
                    continue
                if line.startswith('//'):
                    for i in range(0, len(chem_list)):
                        coefficient_dict[chem_list[i]] = coefficient_list[i]
                    self.dict_rxn_coefficients[reaction] = coefficient_dict
                    coefficient_dict = {}
                    reaction = None
                    chem_list = []
                    coefficient_list = []

        self.check_cofactor_recycling()
        linear_reaction_coefficients(self.model)
        self.model.optimize()
        self.growth_rate = self.model.reactions.get_by_id(self.biomass_rxn).flux
        self.model.reactions.get_by_id(self.biomass_rxn).upper_bound = self.growth_rate * 0.8
        self.model.optimize()
        self.cofactor_recyc_val.append(self.model.reactions.get_by_id(self.ATPS4r).flux)
        self.cofactor_recyc_val.append(self.model.reactions.get_by_id(self.THD2).flux)
        self.cofactor_recyc_val.append(self.model.reactions.get_by_id(self.NADH16).flux)

    def add_rxn(self, list_of_rxns):
        list_of_rxns = list_of_rxns.split("\n")
        self.model_copy = self.model.copy()
        self.model_copy.reactions.get_by_id(self.biomass_rxn).lower_bound = 0.8 * self.growth_rate
        self.dict_metabolite2 = copy.deepcopy(self.dict_metabolite2biocyc)
        isLastRxn = False
        for i in range(len(list_of_rxns)):
            if list_of_rxns[i] == "":
                list_of_rxns.pop(i)
        possible = False

        for i in range(1, len(list_of_rxns) + 1):
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
                    isLastRxn = False
                continue

            self.model_copy.add_reactions({reaction})

            if isLastRxn:
                self.model_copy.objective = reaction.id

            metabolites = tabs[1].split(" --> ")
            reactants = metabolites[0].split(' ')
            products = metabolites[1].split(' ')
            coefficient_dict = self.dict_rxn_coefficients.get(tabs[0])
            for c in reactants:
                if c == '':
                    continue
                if c == "ETF-Oxidized" or c == "ETF-Reduced":
                    continue
                coefficient = coefficient_dict.get(c)
                if "META:" + c in self.dict_metabolite2:
                    c = self.dict_metabolite2.get("META:" + c)
                    reaction.add_metabolites({c.id: -1 * coefficient})
                    continue

                if self.dict_biocyc2metanetx.get(c) in self.dict_metabolite2metanetx:
                    c = self.dict_metabolite2metanetx.get(self.dict_biocyc2metanetx.get(c))
                    reaction.add_metabolites({c.id: -1 * coefficient})
                    continue

                if "ETR-Quinones" in c or "ETR-Quinols" in c or "Acceptor" in c or "Donor-H2" in c or "NADH-P-OR-NOP" in\
                        c or "NAD-P-OR-NOP" in c:
                    print("detected generic compound name. Provide a native " + c)
                    quinone = input()
                    try:
                        self.model_copy.metabolites.get_by_id(quinone)
                    except:
                        print('metabolite not found in model')
                    reaction.add_metabolites({quinone: -1 * coefficient})
                    self.dict_metabolite2biocyc["META:" + c] = self.model.metabolites.get_by_id(quinone)
                    self.dict_metabolite2["META:" + c] = self.model.metabolites.get_by_id(quinone)
                    continue

                else:
                    metabolite = self.add_metabolite("META:" + c)
                    reaction.add_metabolites({metabolite.id: -1 * coefficient})

            for c in products:
                if c == '':
                    continue
                if c == "ETF-Oxidized" or c == "ETF-Reduced":
                    continue
                coefficient = coefficient_dict.get(c)
                if "META:" + c in self.dict_metabolite2:
                    c = self.dict_metabolite2.get("META:" + c)
                    reaction.add_metabolites({c.id: coefficient})
                    continue
                if self.dict_biocyc2metanetx.get(c) in self.dict_metabolite2metanetx:
                    c = self.dict_metabolite2metanetx.get(self.dict_biocyc2metanetx.get(c))
                    reaction.add_metabolites({c.id: coefficient})
                    continue

                if "ETR-Quinones" in c or "ETR-Quinols" in c or "Acceptor" in c or "Donor-H2" in c or "NADH-P-OR-NOP" in c or "NAD-P-OR-NOP" in c:
                    print("detected generic compound name. Provide a native " + c)
                    quinone = input()
                    try:
                        self.model_copy.metabolites.get_by_id(quinone)
                    except:
                        print('metabolite not found in model')
                    reaction.add_metabolites({quinone: coefficient})
                    self.dict_metabolite2biocyc["META:" + c] = self.model.metabolites.get_by_id(quinone)
                    self.dict_metabolite2["META:" + c] = self.model.metabolites.get_by_id(quinone)
                    continue

                else:
                    metabolite = self.add_metabolite("META:" + c)
                    reaction.add_metabolites({metabolite.id: coefficient})
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
            if not self.dict_biocyc2metanetx.get(
                    c) in self.dict_metabolite2metanetx and not "META:" + c in self.dict_metabolite2:
                return False
        return True

    def add_metabolite(self, biocyc_unique_id):
        metabolite = Metabolite(biocyc_unique_id, compartment='c')
        self.dict_metabolite2[biocyc_unique_id] = metabolite
        self.model_copy.add_metabolites(metabolite)
        return metabolite

    def check_cofactor_recycling(self):
        if self.ATPS4r is None:
            reaction = Reaction('ATPS4r')
            self.model.add_reactions([reaction])
            reaction.add_metabolites({'adp_c': -1.0, 'pi_c': -1.0, 'h_e': -4.0, 'atp_c': 1.0, 'h_c': 3.0, 'h2o_c': 1.0})
        if self.THD2 is None:
            reaction = Reaction('THD2')
            self.model.add_reactions([reaction])
            reaction.add_metabolites(
                {'nadh_c': -1.0, 'nadp_c': -1.0, 'h_e': -2.0, 'nad_c': 1.0, 'h_c': 2.0, 'nadph_c': 1.0})

    def run(self, list_of_pathways):
        idx = 0
        data_frame_out = []
        for pathway in list_of_pathways:
            try:
                if not self.add_rxn(pathway):
                    idx = idx + 1
                    continue
            except:
                return data_frame_out

            # Theoretical yield calculation
            linear_reaction_coefficients(self.model)
            solution = self.model_copy.optimize()
            pfba_solution = cobra.flux_analysis.pfba(self.model_copy)
            theoretical_yield = solution.objective_value
            if theoretical_yield == 0.0:
                idx = idx + 1
                continue

            eng_ATPS4r = self.model_copy.reactions.get_by_id(self.ATPS4r).flux
            eng_THD2 = self.model_copy.reactions.get_by_id(self.THD2).flux
            eng_NADH16 = self.model_copy.reactions.get_by_id(self.NADH16).flux

            # cofactor calculation
            # nadh_summary = self.model_copy.metabolites.nadh_c.summary()
            # nadh_consumption = nadh_summary.consuming_flux["flux"].sum()

            # fmnh2_summary = self.model_copy.metabolites.fmnh2_c.summary()
            # fmnh2_consumption = fmnh2_summary.consuming_flux["flux"].sum()

            # fadh2_summary = self.model_copy.metabolites.fadh2_c.summary()
            # fadh2_consumption = fadh2_summary.consuming_flux["flux"].sum()

            # nadph_summary = self.model_copy.metabolites.nadph_c.summary()
            # nadph_consumption = nadph_summary.consuming_flux["flux"].sum()
            #
            # atp_summary = self.model_copy.metabolites.atp_c.summary()
            # atp_consumption = atp_summary.consuming_flux["flux"].sum()

            # FVA span calculation
            FVA = flux_variability_analysis(self.model_copy)
            dif = FVA["maximum"] - FVA["minimum"]
            fva_dif = dif.sum()
            # deletion_results = single_reaction_deletion(self.model_copy)
            # i = 0
            # deletion_set = set()
            # while i <= 3:
            #     index = deletion_results["growth"].idxmax()
            #     max_gene = deletion_results.loc[index, "ids"]
            #     max_gene1 = "".join(max_gene)
            #     if "BIOMASS" in max_gene1 or "ATPM" in max_gene1:
            #         deletion_results.drop(index, axis=0, inplace=True)
            #         continue
            #     max_gene1 = self.model_copy.reactions.get_by_id(max_gene1)
            #     deletion_set.add(max_gene1)
            #     deletion_results.drop(index, axis=0, inplace=True)
            #     i = i + 1
            # #
            # self.model_copy.remove_reactions(deletion_set)
            eng_yield = 0

            anaerobic_medium = self.model.medium
            anaerobic_medium['EX_o2_e'] = 0.0
            self.model_copy.medium = anaerobic_medium
            yield_anaerobic = self.model_copy.slim_optimize()
            fva_dif_anaerobic = 'NaN'
            if type(yield_anaerobic) is not float:
                FVA = flux_variability_analysis(self.model_copy)
                dif = FVA["maximum"] - FVA["minimum"]
                fva_dif_anaerobic = dif.sum()
                self.model_copy.medium = self.model.medium

            entry = [idx, theoretical_yield, yield_anaerobic, self.cofactor_recyc_val[0] - eng_ATPS4r,
                     self.cofactor_recyc_val[1] - eng_THD2, self.cofactor_recyc_val[2] - eng_NADH16,
                     fva_dif, fva_dif_anaerobic, eng_yield]
            data_frame_out.append(entry)
            #
            print('Theoretical_yield: ' + str(theoretical_yield) + '\n' +
                  "yield_anaerobic: " + str(yield_anaerobic) + '\n' +
                  "eng_ATPS4r" + str(eng_ATPS4r) + '\n' +
                  'eng_THD2' + str(eng_THD2) + '\n' +
                  'eng_NAHD16' + str(eng_NADH16) + '\n' +
                  'eng_yield' + str(eng_yield) + '\n' +
                  'FVA_dif' + str(fva_dif) + '\n' +
                  'FVA_dif_anaerobic ' + str(fva_dif_anaerobic) + '\n' +
                  'idx' + str(idx) + '\n' + pathway)

            self.model_copy = None
            idx = idx + 1
        return data_frame_out


def runner(model_path, list_of_paths):
    # model_path = sys.argv[1]
    # list_of_paths = sys.argv[2]
    list_of_paths = list_of_paths.split("//")
    converter = CobraConverter(model_path)
    out = converter.run(list_of_paths)
    display(out)


#

if __name__ == "__main__":
    runner(
        "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/iML1515.xml",
        'RXN-161\tBUTANAL NADH PROTON  --> BUTANOL NAD\n'
        'BUTANAL-DEHYDROGENASE-RXN\tNADH-P-OR-NOP BUTYRYL-COA PROTON  --> CO-A NAD-P-OR-NOP BUTANAL\n'
        'RXN-12558\tCROTONYL-COA NADH PROTON  --> NAD BUTYRYL-COA \n'
        'GLUTARYL-COA-DEHYDROGENASE-RXN\tGLUTARYL-COA ETF-Oxidized PROTON  --> CROTONYL-COA CARBON-DIOXIDE ETF-Reduced\n'
        '2-KETO-ADIPATE-DEHYDROG-RXN\t2K-ADIPATE CO-A NAD  --> CARBON-DIOXIDE GLUTARYL-COA NADH\n'
        '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tCPD-468 2-KETOGLUTARATE  --> 2K-ADIPATE GLT \n'
        '//ENZRXN-201-RXN\tNADPH BUTANAL PROTON  --> NADP BUTANOL \n'
        'RXN0-6973\tOXYGEN-MOLECULE CPD-3744 FMNH2  --> FMN SO3 WATER BUTANAL PROTON \n'
        '//ENZRXN-201-RXN\tNADPH BUTANAL PROTON  --> NADP BUTANOL \n'
        'BUTANAL-DEHYDROGENASE-RXN\tNADH-P-OR-NOP BUTYRYL-COA PROTON  --> CO-A NAD-P-OR-NOP BUTANAL \n'
        'ISOBUTYRYL-COA-MUTASE-RXN\tISOBUTYRYL-COA  --> BUTYRYL-COA \n'
        '1.2.1.25-RXN\tCO-A 2-KETO-ISOVALERATE NAD  --> ISOBUTYRYL-COA CARBON-DIOXIDE NADH \n'
        'VALINE-PYRUVATE-AMINOTRANSFER-RXN\tVAL PYRUVATE  --> L-ALPHA-ALANINE 2-KETO-ISOVALERATE  \n'
    )
