import os
import unittest
import cobra
import pandas as pd
import sbol3
from cobra.flux_analysis import flux_variability_analysis

from CobraConverter import CobraConverter
from Config import Config
from cobra.util.solver import linear_reaction_coefficients

from SBOLDocumenter import SBOLDocumenter


class CobraConverterTest(unittest.TestCase):
    cfg = Config('args.yml')
    converter = CobraConverter('bigg_models/iML1515.xml', cfg)

    def test_stoichiometry(self):
        rxn_str = '1.10.3.7-RXN\tOXYGEN-MOLECULE PROTON --> CPD-1901 WATER\n'
        self.converter.add_rxn(rxn_str)
        rxn = self.converter.model_copy.reactions.get_by_id('1.10.3.7-RXN')
        self.assertEqual(rxn.metabolites[self.converter.model_copy.metabolites.o2_c], -1)
        self.assertEqual(rxn.metabolites[self.converter.model_copy.metabolites.h_c], -2)
        self.assertEqual(rxn.metabolites[self.converter.model_copy.metabolites.h2o_c], 2)
        self.assertEqual(rxn.metabolites[self.converter.model_copy.metabolites.get_by_id('META:CPD-1901')], 2)

    def test_generic_compound(self):
        rxn_str = 'RXN-16062\t1-KETO-2-METHYLVALERATE NAD-P-OR-NOP --> 2-ACETO-2-HYDROXY-BUTYRATE NADH-P-OR-NOP PROTON'
        self.converter.add_rxn(rxn_str)
        rxn = self.converter.model_copy.reactions.get_by_id('RXN-16062')
        self.assertEqual(rxn.metabolites[self.converter.model_copy.metabolites.nadh_c], 1)
        with self.assertRaises(KeyError):
            rxn.metabolites[self.converter.model_copy.metabolites.nadph_c]


    def test_integrated_model(self):
        # Create a new model manually first and compare results to
        model = cobra.io.read_sbml_model('bigg_models/iML1515.xml')
        linear_reaction_coefficients(model)
        solution = model.optimize()
        growth_rate = solution.objective_value
        model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').lower_bound = 0.8 * growth_rate

        propanediol_dehydrogenase = cobra.Reaction('propanediol_dehydrogenase')
        model.add_reactions({propanediol_dehydrogenase})
        propanediol_dehydrogenase.add_metabolites({'h_c': -1, 'nadh_c': -1, '3hppnl_c': -1, '13ppd_c': 1, 'nad_c': 1})

        glycerol_dehydratase = cobra.Reaction('glycerol_dehydratase')
        model.add_reactions({glycerol_dehydratase})
        glycerol_dehydratase.add_metabolites({'glyc_c': -1, '3hppnl_c': 1, 'h2o_c': 1})
        model.add_boundary(model.metabolites.get_by_id('13ppd_c'),
                           type="demand")
        model.objective = propanediol_dehydrogenase
        theoretical_yield_manual = model.optimize().objective_value

        FVA = flux_variability_analysis(model)
        dif = FVA["maximum"] - FVA["minimum"]
        fva_dif_manual = dif.sum()

        self.converter.add_rxn('13-PROPANEDIOL-DEHYDROGENASE-RXN_rev\tPROTON HYDROXYPROPANAL NADH  --> CPD-347 NAD\n' \
                               'GLYCEROL-DEHYDRATASE-RXN\tGLYCEROL  --> HYDROXYPROPANAL WATER\n')

        solution = self.converter.model_copy.optimize().objective_value
        FVA = flux_variability_analysis(self.converter.model_copy)
        dif = FVA["maximum"] - FVA["minimum"]
        fva_dif = dif.sum()
        self.assertTrue(abs(fva_dif - fva_dif_manual) // fva_dif < 0.0001)
        self.assertTrue(abs(solution - theoretical_yield_manual) // solution < 0.0001)

    def test_sbol_validation(self):
        # validate completeness of SBOL files.a SBOL document is self-contained or 'complete' if every SBOL object
        # referred to in teh document is contained in the document.
        dirname = os.path.dirname(__file__)
        rxn_dat_path = os.path.join(dirname, 'data', 'reactions.txt')
        chem_dat_path = os.path.join(dirname, 'data', 'chems.txt')
        file_writer = SBOLDocumenter(rxn_dat_path, chem_dat_path, self.cfg.get('save_SBOL_files_to'))

        out = self.converter.run(['13-PROPANEDIOL-DEHYDROGENASE-RXN_rev\tPROTON HYDROXYPROPANAL NADH  --> CPD-347 NAD\n' \
                               'GLYCEROL-DEHYDRATASE-RXN\tGLYCEROL  --> HYDROXYPROPANAL WATER\n'])
        df = pd.DataFrame(out,
                      columns=['idx', 'theoretical_yield', 'eng_atp', 'eng_nad', 'eng_nadp', 'fva_dif',
                               'yield_anaerobic', 'anaerobic_atp_use', 'anaerobic_nadh_use', 'anaerobic_nadph_use',
                               'fva_dif_anaerobic', 'model'])
        for row in df.iterrows():
            idx = row[1].T.idx
            path = '13-PROPANEDIOL-DEHYDROGENASE-RXN_rev\tPROTON HYDROXYPROPANAL NADH  --> CPD-347 NAD\n' \
                               'GLYCEROL-DEHYDRATASE-RXN\tGLYCEROL  --> HYDROXYPROPANAL WATER\n'
            file_writer.add_new_path(path, row, idx)
        result_path = os.path.join(dirname, 'results')
        for filename in os.listdir(result_path):
            f = os.path.join(result_path, filename)
            doc = sbol3.Document()
            doc.read(f)
            report = doc.validate()
            self.assertEqual(len(report), 0)

if __name__ == '__main__':
    unittest.main()
