import unittest
import subprocess

from IPython.core.display import display

from CobraConverter import CobraConverter

from main import ranker
import pandas as pd


class Test(unittest.TestCase):

    def test(self):
        """
        Test that it can sum a list of integers
        """
        print('Enter target id, precursor id (optional), and max path length separated by space: ')
        out = subprocess.Popen(["java", "-jar",
                                "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/PathEnumerator.jar"],
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output_itr = iter(out.stdout.readline, b'')
        out_lst = []
        for line in output_itr:
            line = line.decode('utf-8')
            out_lst.append(line)
        out_lst.pop(0)
        list_of_paths = ''.join(out_lst)
        list_of_paths = list_of_paths.split("//")
        print("Enter path to chassis model")
        model_path = input()

        converter = CobraConverter(model_path)
        print('for carbon fixation? ')
        carbon_fixation = input()

        if carbon_fixation == "carbon_fixation":
            medium = converter.model.medium
            medium['EX_photon_e'] = 100
            # medium['EX_glc__D_e'] = 0
            converter.model.reactions.EX_photon_e.lower_bound = -100
            converter.model.reactions.EX_glc__D_e.lower_bound = 0
            converter.model.reactions.EX_co2_e.lower_bound = -3.7
            # converter.model.reactions.EX_hco3_e.lower_bound = -3.7
            converter.model.medium = medium

        df_output = converter.run(list_of_paths)
        df = pd.DataFrame(df_output,
                          columns=['index', 'Theoretical_yield', 'yield_anaerobic', 'eng_ATPM', 'eng_ATPS4r',
                                   'eng_NADH16', 'fva_dif', 'fva_dif_anaerobic', 'eng_yield'])
        rankingparameter = input("Enter ranking parameter: " + "\n"
                                 "0: idx" + "\n"
                                 "1: descending order of theoretical yield" + "\n"
                                 "2: ascending order of FVA" + "\n"
                                 "3: descending order of ATP use" + "\n"
                                 "4: descending order of NADH use" + "\n"
                                 "5: descending order of NADPH use")
        rankingparameter = int(rankingparameter)
        ranked = ranker(df, rankingparameter)
        display(ranked)
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
