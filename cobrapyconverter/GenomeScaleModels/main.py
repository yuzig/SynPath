import os
import shutil
import sys
import subprocess
from flask import Flask, render_template, request, send_file

from IPython.core.display import display
from werkzeug.utils import secure_filename

from GenomeScaleModels import Escher_visualizer
from cobraConverterFromFile import cobraConverterFromFile
import pandas as pd

from Config import Config
from SBOLDocumenter import SBOLDocumenter
import d3flux as d3f
from cobra.core import Model

app = Flask(__name__)

list_of_paths = []
filename = None
@app.route('/')
def form():
    return render_template('form.html')

@app.route('/intermediate_step', methods=['POST', 'GET'])
def intermediate_step():
    global filename
    if request.method == 'POST':
        target = request.form['target_id']
        precursor = request.form['precursor_id']
        max_len = request.form['max_path_length']
        file = request.files['myfile']

        if file:
            filename = secure_filename(file.filename)
            file.save(filename)

        dirname = os.path.dirname(__file__)

        lst_input = [target, precursor, max_len]
        lst_input = [i for i in lst_input if i]

        PathEnumerator_jar = os.path.join(dirname, 'PathEnumerator.jar')
        out = subprocess.Popen(["java", "-jar", PathEnumerator_jar] + lst_input,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output_itr = iter(out.stdout.readline, b'')

        out_lst = []
        for line in output_itr:
            line = line.decode('utf-8')
            out_lst.append(line)
        out.stdout.close()
        global list_of_paths
        list_of_paths = []
        try:
            i = out_lst.index("//\n")
            del out_lst[i]
        except:
            pass
        list_of_paths = ''.join(out_lst)
        list_of_paths = list_of_paths.split("//")
        total_n_output = len(list_of_paths)
        return render_template('intermediate_step.html',msg=total_n_output)

@app.route('/data', methods=['POST', 'GET'])
def data():
    if request.method == 'GET':
        return f"The URL /data is accessed directly. Try going to '/form' to submit form"
    if request.method == 'POST':
        dirname = os.path.dirname(__file__)
        cfg_path = os.path.join(dirname, 'args.yml')
        cfg = Config(cfg_path)

        converter = cobraConverterFromFile(filename, cfg)
        df_output = converter.run(list_of_paths)
        df = pd.DataFrame(df_output,
                          columns=['idx', 'theoretical_yield', 'eng_atp', 'eng_nad', 'eng_nadp', 'fva_dif','o2_use',
                                   'yield_anaerobic', 'anaerobic_atp_use', 'anaerobic_nadh_use', 'anaerobic_nadph_use',
                                   'fva_dif_anaerobic', 'model'])
        rankingparameter = cfg.get('rank_param')
        ranked = ranker(df, rankingparameter)
        lst_paths = converter.get_pathways()
        os.remove(filename)
        rxn_dat_path = os.path.join(dirname, 'data','reactions.txt')
        chem_dat_path = os.path.join(dirname, 'data','chems.txt')
        file_writer = SBOLDocumenter(rxn_dat_path, chem_dat_path, "results")

        map = converter.get_all_added_rxns()
        escher_map = visualizer(map)

        for row in df.iterrows():
            idx = row[1].T.idx
            paths = list_of_paths[idx]
            file_writer.add_new_path(paths, row, idx)

        if os.path.exists('results.zip'):
            os.remove('results.zip')

        here = os.path.dirname(os.path.realpath(__file__))
        result_directory = os.path.join(here, "results")
        shutil.make_archive('results', 'zip', result_directory)

        return render_template('data.html', lst_paths=lst_paths, column_names=ranked.columns.values,
                               row_data=list(ranked.values.tolist()), zip=zip, map=escher_map)
def visualizer(reactions_set):

    model = Model("escher_map")
    for reaction in reactions_set:
        try:
            model.reactions.get_by_id(reaction.id)
        except:
            model.add_reactions([reaction])
    list_of_common_cofactors = ['nadh_c','nad_c','nadp_c','nadph_c','h_c', 'atp_c','adp_c',
                                 'amp_c','coa_c','h2o','o2_c','h2_c','co2_c','q8_c','q8h2_c',
                                 'pi_c', 'nh4_c']
    for c in list_of_common_cofactors:
        try:
            d3f.update_cofactors(model, [c])
        except:
            continue
    return d3f.flux_map(model)

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
#
@app.route('/download')
def download_file():
    path = "results.zip"
    return send_file(path, as_attachment=True)

@app.route('/download_chem_data')
def download_chem_data():
    path = "data/chems.txt"
    return send_file(path, as_attachment=True)

if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=5000)

