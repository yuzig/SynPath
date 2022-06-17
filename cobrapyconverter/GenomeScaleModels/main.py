import subprocess

from IPython.core.display import display
from CobraConverter import CobraConverter
import pandas as pd

from SBOLDocumenter import SBOLDocumenter


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
        return df.sort_values(by=['theoretical_yield'], ascending=False)
    if order == 2:
        return df.sort_values(by=['fva_dif'])
    if order == 3:
        return df.sort_values(by=['yield_anaerobic'])
    if order == 4:
        return df.sort_values(by=['fva_dif_anaerobic'])


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
                      columns=['idx', 'theoretical_yield', 'eng_atp', 'eng_nad', 'eng_nadp', 'fva_dif',
                               'yield_anaerobic', 'anaerobic_atp_use', 'anaerobic_nadh_use', 'anaerobic_nadph_use',
                               'fva_dif_anaerobic', 'model'])
    print("Enter ranking parameter: " + "\n"
          "0: idx" + "\n"
          "1: descending order of theoretical yield" + "\n"
          "2: ascending order of FVA_aerobic" + "\n"
          "3: descending order of anaerobic theoretical yield" + "\n"
          "4: descending order of anaerobic FVA span" + "\n")
    rankingparameter = input()
    rankingparameter = int(rankingparameter)
    ranked = ranker(df, rankingparameter)
    display(ranked)

    print('Convert results to SBOL files?')
    convert_to_SBOLf = input()
    if convert_to_SBOLf == "yes":
        rxn_dat_path = '/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/reactions.txt'
        chem_dat_path = '/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/chems.txt'
        file_writer = SBOLDocumenter(rxn_dat_path, chem_dat_path)
        for params in ranked:
            idx = params[0]
            path = list_of_paths[idx]
            file_writer.add_new_path(path, params)
