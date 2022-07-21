import shutil

import sbol3
import os

from data_models.Chemicals import Chemicals


class SBOLDocumenter:
    """
    SBOLDocumenter converts pathways and analysis results into SBOL compliant format
    Each SBOLdocument is comprised of only one pathway
    """
    dict_reaction = None
    dict_chemicals = None
    result_directory = None


    def __init__(self, rxn_dat_file, chem_dat_file, result_directory):
        dict_chemicals = {}
        dict_reaction = {}

        with open(chem_dat_file) as f:
            lines = f.readlines()
            inchi = ''
            unique_id = ''
            name = ''
            smiles = ''
            for line in lines:
                if line.startswith('UNIQUE-ID'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    unique_id = line[1]
                    continue
                if line.startswith('COMMON-NAME'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    name = line[1]
                    continue
                if line.startswith('INCHI '):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    inchi = line[1]
                    continue
                if line.startswith('SMILES'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    smiles = line[1]
                    continue
                if line.startswith('//'):
                    compound = Chemicals(unique_id, name, inchi, smiles)
                    dict_chemicals[unique_id] = compound
                    inchi = ''
                    unique_id = ''
                    name = ''
                    smiles = ''
            self.dict_chemicals = dict_chemicals

        with open(rxn_dat_file) as f:
            lines = f.readlines()
            ecnum = ''
            unique_id = ''
            for line in lines:
                if line.startswith('UNIQUE-ID'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    unique_id = line[1]
                    continue
                if line.startswith('EC-NUMBER'):
                    line = line.replace("\n", "")
                    line = line.split(" - ")
                    ecnum = line[1]
                    dict_reaction[unique_id] = ecnum
                    continue
                if line.startswith('//'):
                    unique_id = ''
                    ecnum = ''
                    continue
            self.dict_reaction = dict_reaction

        here = os.path.dirname(os.path.realpath(__file__))
        self.result_directory = os.path.join(here, result_directory)
        if os.path.exists(self.result_directory):
            shutil.rmtree(self.result_directory)
        os.mkdir(os.path.join(here, self.result_directory))

    def add_new_path(self, reactions_str, pathway_param, idx):

        lines = reactions_str.split('\n')

        namespace = 'http://sbolstandard.org/testfiles' + str(idx)
        sbol3.set_namespace(namespace)
        doc = sbol3.Document()
        intro = sbol3.Collection('pathway')
        description = 'theoretical_yield = ' + str(pathway_param[1].T.theoretical_yield) + '//' + "eng_atp = " + str(pathway_param[1].T.eng_atp) + "//" + 'nad_aerobic = ' + str(pathway_param[1].T.eng_nad) + '//' + 'nadp_aerobic = ' + str(pathway_param[1].T.eng_nadp) + '//' + 'fva_dif_aerobic = ' + str(pathway_param[1].T.fva_dif) + '//' + 'yield_anaerobic = ' + str(pathway_param[1].T.yield_anaerobic) + '//' + 'anaerobic_atp_use = ' + str(pathway_param[1].T.anaerobic_atp_use) + '//' + 'anaerobic_nad_use = ' + str(pathway_param[1].T.anaerobic_nadh_use) + '//' + 'anaerobic_nadp_use = ' + str(pathway_param[1].T.anaerobic_nadph_use) + '//' + 'fva_dif_anaerobic = ' + str(pathway_param[1].T.fva_dif_anaerobic ) + '//' + 'index = ' + str(pathway_param[0]) + '//' + 'model = ' + pathway_param[1].T.model
        intro.description = description
        doc.add(intro)
        i = 0
        for line in lines:
            if line == '':
                continue
            self.add_rxn(line, doc, namespace, i, intro)
            i += 1
        fileName = 'pathway' + str(idx) + '.xml'

        filepath = os.path.join(self.result_directory,fileName)
        doc.write(filepath)
        print('pathway ' + str(idx) + ' successfully converted to SBOL file to ' + self.result_directory)

    def add_rxn(self, rxn_str, SBOL_doc, namespace, i, intro):
        line = rxn_str.split('\t')
        reaction = line[0]
        reaction= reaction.replace('-', '_')
        reaction = reaction.replace('.', '_')
        reaction = reaction.replace('+','')
        reaction_str = 'r_' + reaction
        reaction_component = sbol3.Component(reaction_str, sbol3.SBO_BIOCHEMICAL_REACTION)
        reaction_component.description = self.dict_reaction.get(reaction)
        intro.members += [reaction_component]
        SBOL_doc.add(reaction_component)
        interaction = sbol3.Interaction('rxn' + str(i))
        interaction.type = sbol3.SBO_BIOCHEMICAL_REACTION
        chemicals = line[1].split('-->')
        reactants = chemicals[0].split(' ')
        products = chemicals[1].split(' ')

        if type(reactants) == str:
            reactants = [reactants]

        if type(products) == str:
            products = [products]

        for reactant in reactants:
            if reactant == '':
                continue
            reactant_str = reactant.replace(' ', '')
            reactant_str = reactant_str.replace('-', '_')
            if SBOL_doc.find(namespace + '/' + 'c_' + reactant_str) is not None:
                reactant_component = SBOL_doc.find(namespace + '/' + 'c_' + reactant_str)
            else:
                reactant_component = sbol3.Component('c_' + reactant_str, sbol3.SBO_SIMPLE_CHEMICAL)
                compound = self.dict_chemicals.get(reactant)
                if compound:
                    seq = sbol3.Sequence('seq_'+reactant_str, elements=compound.inchi, encoding=sbol3.INCHI_ENCODING)
                    reactant_component.sequences = [seq]
                if SBOL_doc.find(namespace + '/' + 'seq_' + reactant_str) is None:
                    try:
                        SBOL_doc.add(seq)
                    except ValueError:
                        pass

            reactant_sc = sbol3.SubComponent(reactant_component, roles=sbol3.SBO_REACTANT)
            reaction_component.features += [reactant_sc]

        for product in products:
            if product == '':
                continue
            product = product.replace(' ', '')
            product_str = product.replace('-', '_')
            if SBOL_doc.find(namespace + '/' + 'c_' + product_str) is not None:
                product_component = SBOL_doc.find(namespace + '/' + 'c_' + product_str)
            else:
                product_component = sbol3.Component('c_' + product_str, sbol3.SBO_SIMPLE_CHEMICAL)
                compound = self.dict_chemicals.get(product)
                if compound:
                    seq = sbol3.Sequence('seq_' + product_str, elements=compound.inchi, encoding=sbol3.INCHI_ENCODING)
                    product_component.sequences = [seq]
                if SBOL_doc.find(namespace + '/' + 'seq_' + product_str) is None:
                    try:
                        SBOL_doc.add(seq)
                    except ValueError:
                        pass

            product_sc = sbol3.SubComponent(product_component, roles=sbol3.SBO_PRODUCT)
            reaction_component.features += [product_sc]

        reaction_component.interactions = [interaction]

