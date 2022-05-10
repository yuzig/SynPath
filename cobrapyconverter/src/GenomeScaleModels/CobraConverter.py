import cobra
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
import copy


class CobraConverter:
    dict_metabolite2biocyc = {}
    dict_metabolite2metanetx = {}
    dict_biocyc2metanetx = {}
    dict_rxn = {}
    dict_metabolite2 = {}
    model = None
    model_copy = None

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

        with open('biocyc2metanetx.txt') as f:
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
            if r == '':
                continue
            tabs = r.split("\t")

            if "META:" + tabs[0] in self.dict_rxn or tabs[0] == '':
                if isLastRxn:
                    self.model_copy.objective = self.dict_rxn.get("META:" + tabs[0]).id
                    isLastRxn = False
                continue

            reaction = Reaction(tabs[0])
            self.model_copy.add_reactions({reaction})

            if isLastRxn:
                self.model_copy.objective = reaction.id

            metabolites = tabs[1].split(" --> ")
            reactants = metabolites[0].split(' ')
            products = metabolites[1].split(' ')
            for c in reactants:
                if c == '':
                    continue
                if "META:" + c in self.dict_metabolite2:
                    c = self.dict_metabolite2.get("META:" + c)
                    reaction.add_metabolites({c.id: -1})
                    continue;

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
        theoretical_yield = []
        for pathway in list_of_pathways:
            self.add_rxn(pathway)
            linear_reaction_coefficients(self.model)
            self.model_copy.optimize()
            RXN14985 = self.model_copy.reactions.get_by_id('RXN-14985')
            ETHYLMALATESYNTHASERXN = self.model_copy.reactions.get_by_id('3-ETHYLMALATE-SYNTHASE-RXN')
            RXN12558 = self.model_copy.reactions.get_by_id('RXN-12558')
            ENZRX201RXN = self.model_copy.reactions.get_by_id('ENZRXN-201-RXN')
            RXN18210 = self.model_copy.reactions.get_by_id('RXN-18210')
            RXN18211 = self.model_copy.reactions.get_by_id('RXN-18211')
            glucose_exchange = self.model_copy.reactions.get_by_id('EX_glc__D_e')

            rxns = [ENZRX201RXN, ETHYLMALATESYNTHASERXN, RXN14985,RXN18210,RXN18211, RXN12558,glucose_exchange]
            out = flux_variability_analysis(self.model_copy, rxns)

            print(out)
            # linear_reaction_coefficients(self.model)
            # model_solution = self.model_copy.optimize()
            # theoretical_yield.append(model_solution.objective_value)
        #     self.model_copy = None
        #     self.dict_metabolite2 = {}
        # print(theoretical_yield)
        # return theoretical_yield


class runCobrapy:
    def runner(model_path, list_of_paths):
        list_of_paths = list_of_paths.split('//')
        converter = CobraConverter(model_path)
        converter.run(list_of_paths)


if __name__ == "__main__":
    runCobrapy.runner(
        "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/src/GenomeScaleModels/iJO1366.xml",
       "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n"
       "RXN-14985\tCPD-3618 PROTON  --> CARBON-DIOXIDE BUTANAL\n"
       "RXN-18211\tCPD-19492 PROTON  --> CPD-3618 CARBON-DIOXIDE\n"
       "RXN-18210\tNAD CPD-1130  --> NADH CPD-19492 PROTON\n"
       "3-ETHYLMALATE-SYNTHASE-RXN\tWATER GLYOX BUTYRYL-COA  --> CO-A CPD-1130 PROTON\n"
       "RXN-12558\tNADH CROTONYL-COA PROTON  --> NAD BUTYRYL-COA\n"
    )

# columns= [["ID", "name", "Biocyc_Annotation"]]
# for m in model.metabolites:
#     chem = [m.id, m.name, m.annotation.get('biocyc')]
#     chem = m.annotation.get('biocyc')
#
# columns.append(chem)
# df = pd.DataFrame(columns)
# df.to_csv('ijo1366_metabolites.csv')
#
# def run(string):


# #
# columns_rxn= [["ID", "name", "ec_num","Biocyc_Annotation"]]
# for r in model.reactions:
#     if r.annotation.get('biocyc') is None:
#         continue
#     rxn = [r.id, r.name, r.annotation.get('ec-code'), r.annotation.get('biocyc')]
#     columns_rxn.append(rxn)
# df = pd.DataFrame(columns_rxn)
# df.to_csv('ijo1366_reactions.csv')
#
#
#
#
#
