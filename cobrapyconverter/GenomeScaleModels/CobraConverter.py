import cobra
from IPython.core.display import display
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis
import copy
import math


class CobraConverter:
    dict_metabolite2biocyc = {}
    dict_metabolite2metanetx = {}
    dict_biocyc2metanetx = {}
    dict_rxn = {}
    dict_metabolite2 = {}
    dict_rxn_coefficients = {}
    model = None
    model_copy = None
    growth_rate = 0
    native_atp = 0
    native_nadp = 0
    native_nad = 0
    biomass_rxn = None

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

        for c in self.model.metabolites:
            if c.annotation.get('biocyc') is not None and type(c.annotation.get('biocyc')) == str:
                self.dict_metabolite2biocyc[c.annotation.get('biocyc')] = c
            if c.annotation.get('biocyc') is not None and type(c.annotation.get('biocyc')) == list:
                for biocyc_id in c.annotation.get('biocyc'):
                    self.dict_metabolite2biocyc[biocyc_id] = c
            if c.annotation.get('metanetx.chemical') is not None:
                self.dict_metabolite2metanetx[c.annotation.get('metanetx.chemical')] = c
        with open(
                '/GenomeScaleModels/data/biocyc2metanetx.txt') as f:
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
                '/GenomeScaleModels/data/biocyc2bigg.txt') as f:
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

        with open(
                '/GenomeScaleModels/data/reactions.txt') as f:
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

        linear_reaction_coefficients(self.model)
        self.model.optimize()
        self.growth_rate = self.model.reactions.get_by_id(self.biomass_rxn).flux
        self.model.reactions.get_by_id(self.biomass_rxn).upper_bound = 0.8 * self.growth_rate
        self.model.optimize()
        self.native_atp = self.model.metabolites.atp_c.summary().consuming_flux['flux'].sum()
        self.native_nad = self.model.metabolites.nad_c.summary().consuming_flux['flux'].sum()
        self.native_nadp = self.model.metabolites.nadp_c.summary().consuming_flux['flux'].sum()

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
            tabs[0] = tabs[0].replace('_rev', '')
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

                if "ETR-Quinones" in c or "ETR-Quinols" in c or "Acceptor" in c or "Donor-H2" in c or "NADH-P-OR-NOP" in \
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
            linear_reaction_coefficients(self.model_copy)
            solution = self.model_copy.optimize()
            theoretical_yield = solution.objective_value

            if theoretical_yield == 0.0:
                idx = idx + 1
                continue
            try:
                eng_atp = self.model_copy.metabolites.atp_c.summary().consuming_flux['flux'].sum()
                eng_nad = self.model_copy.metabolites.nad_c.summary().consuming_flux['flux'].sum()
                eng_nadp = self.model_copy.metabolites.nadp_c.summary().consuming_flux['flux'].sum()
            except:
                eng_atp = 'NaN'
                eng_nad = "NaN"
                eng_nadp = "NaN"

            # FVA span calculation
            FVA = flux_variability_analysis(self.model_copy)
            dif = FVA["maximum"] - FVA["minimum"]
            fva_dif = dif.sum()

            anaerobic_medium = self.model.medium
            anaerobic_medium['EX_o2_e'] = 0.0
            self.model_copy.medium = anaerobic_medium
            self.model_copy.reactions.EX_o2_e.lower_bound = 0
            yield_anaerobic = self.model_copy.slim_optimize()
            fva_dif_anaerobic = 'NaN'
            anaerobic_atp_use = 'NaN'
            anaerobic_nadh_use = 'NaN'
            anaerobic_nadph_use = 'NaN'
            if not math.isnan(yield_anaerobic):
                FVA = flux_variability_analysis(self.model_copy)
                dif = FVA["maximum"] - FVA["minimum"]
                fva_dif_anaerobic = dif.sum()
                anaerobic_atp_use = self.model_copy.metabolites.atp_c.summary().consuming_flux['flux'].sum()
                anaerobic_nadh_use = self.model_copy.metabolites.nadh_c.summary().consuming_flux['flux'].sum()
                anaerobic_nadph_use = self.model_copy.metabolites.nadph_c.summary().consuming_flux['flux'].sum()

            entry = [idx, theoretical_yield, eng_atp, eng_nad, eng_nadp, fva_dif,
                     yield_anaerobic, anaerobic_atp_use, anaerobic_nadh_use, anaerobic_nadph_use, fva_dif_anaerobic,
                     self.model.id]
            data_frame_out.append(entry)
            #
            print('Theoretical_yield: ' + str(theoretical_yield) + '\n' +
                  'Aerobic_atp: ' + str(eng_atp) + '\n' +
                  'Aerobic_nadh: ' + str(eng_nad) + '\n' +
                  'Aerobic_nadph: ' + str(eng_nadp) + '\n' +
                  'Aerobic_FVA_span: ' + str(fva_dif) + '\n' +
                  "yield_anaerobic: " + str(yield_anaerobic) + '\n' +
                  'Anaerobic_atp: ' + str(anaerobic_atp_use) + '\n' +
                  'Anaerobic_nadh: ' + str(anaerobic_nadh_use) + '\n' +
                  'Anaerobic_nadph: ' + str(anaerobic_nadph_use) + '\n' +
                  'Anaerobic_FVA_span: ' + str(fva_dif_anaerobic) + '\n' +
                  'idx: ' + str(idx))
            print(pathway)

            self.model_copy = None
            idx = idx + 1
        return data_frame_out


def runner(model_path, list_of_paths):
    list_of_paths = list_of_paths.split("//")
    converter = CobraConverter(model_path)
    out = converter.run(list_of_paths)
    display(out)


