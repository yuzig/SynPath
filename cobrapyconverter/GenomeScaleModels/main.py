import os
import sys
import subprocess

from IPython.core.display import display
from CobraConverter import CobraConverter
import pandas as pd

from Config import Config
from SBOLDocumenter import SBOLDocumenter


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
    dirname = os.path.dirname(__file__)
    cfg_path = os.path.join(dirname, 'args.yml')
    cfg = Config(cfg_path)
    target = cfg.get('target')
    precursor = cfg.get('precursor')
    max_len = str(cfg.get('max_path_len'))
    target = input('target ID: ')
    max_len = input('max len: ')
    precursor = input('precursor: ')
    model_path = input('enter model name: ')
    lst_input = [target, precursor, max_len]
    lst_input = [i for i in lst_input if i]

    dirname = os.path.dirname(__file__)
    PathEnumerator_jar = os.path.join(dirname, 'PathEnumerator.jar')
    out = subprocess.Popen(["java", "-jar", PathEnumerator_jar] + lst_input,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output_itr = iter(out.stdout.readline, b'')

    out_lst = []
    for line in output_itr:
        line = line.decode('utf-8')
        out_lst.append(line)

    out.stdout.close()
    list_of_paths = ''.join(out_lst)
    list_of_paths = list_of_paths.split("//")
    # model_path = cfg.get('bigg_model') + '.xml'
    file_path = os.path.join(dirname, 'bigg_models', model_path)
    converter = CobraConverter(file_path, cfg)

    carbon_fixation = cfg.get('carbon_fixation')

    if carbon_fixation:
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


    rankingparameter = cfg.get('rank_param')
    ranked = ranker(df, rankingparameter)
    display(ranked)

    if cfg.get('to_SBOL_file'):
        rxn_dat_path = os.path.join(dirname, 'data','reactions.txt')
        chem_dat_path = os.path.join(dirname, 'data','chems.txt')
        file_writer = SBOLDocumenter(rxn_dat_path, chem_dat_path, cfg.get('save_SBOL_files_to'))

        for row in df.iterrows():
            idx = row[1].T.idx
            path = list_of_paths[idx]
            file_writer.add_new_path(path, row, idx)
        print('task completed')