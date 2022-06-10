import subprocess

from IPython.core.display import display
from CobraConverter import CobraConverter
import pandas as pd


def pathEnumeratorReader():
    out = subprocess.Popen(["java", "-jar",
                            "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/PathEnumerator.jar"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    subprocess.Popen.kill(out)
    return iter(out.stdout.readline, b'')

def ranker(df, order):
    """
     0: idx
     1: descending order of theoretical yield
     2: ascending order of FVA
     3: descending order of ATP use
     4: descending order of NADH use
     5: descending order of NADPH use
    """
    if order == 0:
        return df
    if order == 1:
        return df.sort_values(by=['Theoretical_yield'], ascending=False)
    if order == 2:
        return df.sort_values(by=['fva_dif'])


if __name__ == "__main__":
    output_itr = pathEnumeratorReader()
    out_lst = []
    for line in output_itr:
        line = line.decode('utf-8')
        out_lst.append(line)
    out_lst.pop(0)
    list_of_paths = ''.join(out_lst)

    list_of_paths = list_of_paths.split("//")
    print("Enter path to chassis models: ")
    model_path = input()

    converter = CobraConverter(model_path)
    print('for carbon fixation? ')
    carbon_fixation = input()

    if carbon_fixation == "yes" or carbon_fixation == 'Y' or carbon_fixation == "y":
        medium = converter.model.medium
        medium['EX_photon_e'] = 100
        converter.model.reactions.EX_photon_e.lower_bound = -100
        converter.model.reactions.EX_glc__D_e.lower_bound = 0
        converter.model.reactions.EX_co2_e.lower_bound = -3.7
        converter.model.reactions.EX_hco3_e.lower_bound = -3.7
        converter.model.medium = medium

    df_output = converter.run(list_of_paths)
    df = pd.DataFrame(df_output,
                      columns=['index', 'Theoretical_yield', 'yield_anaerobic',  'ATPS4r',  'THD2', 'NADH16',
                               'fva_dif', 'fva_dif_anaerobic',
                               'eng_yield'])
    rankingparameter = input("Enter ranking parameter: ")
    rankingparameter = int(rankingparameter)
    ranked = ranker(df, rankingparameter)
    display(ranked)
