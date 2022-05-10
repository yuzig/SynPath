import cobra
from cobra import Model, Reaction, Metabolite
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
        for pathway in list_of_pathways:
            self.add_rxn(pathway)
            linear_reaction_coefficients(self.model)
            self.model_copy.optimize()
            theoretical_yield = flux_variability_analysis(self.model_copy, self.last_reactions[-1])

            print(theoretical_yield)


def runner():
    model_path = sys.argv[1]
    list_of_paths = sys.argv[2]
    list_of_paths = list_of_paths.split("//")
    converter = CobraConverter(model_path)
    converter.run(list_of_paths)


runner()

if __name__ == "__main__":
    runner(
        "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/src/GenomeScaleModels/iJO1366.xml",
        # "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n"
        # "RXN-14985\tCPD-3618 PROTON  --> CARBON-DIOXIDE BUTANAL\n"
        # "RXN-18211\tCPD-19492 PROTON  --> CPD-3618 CARBON-DIOXIDE\n"
        # "RXN-18210\tNAD CPD-1130  --> NADH CPD-19492 PROTON\n"
        # "3-ETHYLMALATE-SYNTHASE-RXN\tWATER GLYOX BUTYRYL-COA  --> CO-A CPD-1130 PROTON\n"
        # "RXN-12558\tNADH CROTONYL-COA PROTON  --> NAD BUTYRYL-COA\n" + "//" +
        # "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n"
        # "RXN-14985\tCPD-3618 PROTON  --> CARBON-DIOXIDE BUTANAL\n "
        # "RXN-14986\tNAD CPD-1130  --> NADH CPD-3618 CARBON-DIOXIDE\n"
        # "3-ETHYLMALATE-SYNTHASE-RXN\tWATER GLYOX BUTYRYL-COA  --> CO-A CPD-1130 PROTON\n"
        # "RXN-12558\tNADH CROTONYL-COA PROTON  --> NAD BUTYRYL-COA\n" + "//" +
        # "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n"
        # "BUTANAL-DEHYDROGENASE-RXN\tBUTYRYL-COA PROTON NADPH  --> CO-A NADP BUTANAL\n" + "//" +
        # "ENZRXN-201-RXN\tNADPH PROTON BUTANAL  --> NADP BUTANOL\n"
        # "RXN-14985\tCPD-3618 PROTON  --> CARBON-DIOXIDE BUTANAL\n"
        # "RXN-14986\tNAD CPD-1130  --> NADH CPD-3618 CARBON-DIOXIDE\n"
        # "3-ETHYLMALATE-SYNTHASE-RXN\tWATER GLYOX BUTYRYL-COA  --> CO-A CPD-1130 PROTON\n"
        # "BUTYRYL-COA-DEHYDROGENASE-RXN\tETF-Reduced CROTONYL-COA  --> ETF-Oxidized BUTYRYL-COA PROTON\n" + "//" +
        # "RXN-12595\tPROTON CPD-13555 NADH  --> CPD-13560 NAD\n"
        # "RXN-12594\tPROTON 4-HYDROXY-BUTYRYL-COA NADH  --> CPD-13555 NAD CO-A\n "
        # "RXN-8889\tACETYL-COA 4-HYDROXY-BUTYRATE  --> ACET 4-HYDROXY-BUTYRYL-COA\n"
        # "4-HYDROXYBUTYRATE-DEHYDROGENASE-RXN\tSUCC-S-ALD PROTON NADH  --> 4-HYDROXY-BUTYRATE NAD\n"
        # "RXN-13328\tGLYOX 4-AMINO-BUTYRATE  --> SUCC-S-ALD GLY\n"
        # "1.5.1.35-RXN\tWATER NAD CPD-6124  --> PROTON 4-AMINO-BUTYRATE NADH\n"
        # "2.6.1.82-RXN\t2-KETOGLUTARATE PUTRESCINE  --> WATER GLT CPD-6124\n"
        # "ORNDECARBOX-RXN\tPROTON L-ORNITHINE  --> CARBON-DIOXIDE PUTRESCINE\n" + "//" +
        "RXN-12595\tPROTON CPD-13555 NADH  --> CPD-13560 NAD\n"
        "RXN-12594\tPROTON 4-HYDROXY-BUTYRYL-COA NADH  --> CPD-13555 NAD CO-A\n"
        "RXN-9092\t4-HYDROXY-BUTYRATE ATP CO-A  --> PPI 4-HYDROXY-BUTYRYL-COA AMP\n"
        "4-HYDROXYBUTYRATE-DEHYDROGENASE-RXN\tSUCC-S-ALD PROTON NADH  --> 4-HYDROXY-BUTYRATE NAD\n"
        "RXN-13328\tGLYOX 4-AMINO-BUTYRATE  --> SUCC-S-ALD GLY\n"
        "GUANIDINOBUTYRASE-RXN\tWATER CPD-592  --> UREA 4-AMINO-BUTYRATE\n"
        "GUANIDINOBUTANAMIDE-NH3-RXN\t4-GUANIDO-BUTYRAMIDE WATER  --> CPD-592 AMMONIUM\n"
        "ARGININE-2-MONOOXYGENASE-RXN\tOXYGEN-MOLECULE ARG  --> 4-GUANIDO-BUTYRAMIDE WATER CARBON-DIOXIDE\n"
    )
