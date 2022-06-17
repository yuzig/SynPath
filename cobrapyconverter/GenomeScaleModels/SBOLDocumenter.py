import sbol3

from Chemicals import Chemicals


class SBOLDocumenter:
    """
    SBOLDocumenter converts pathways and analysis results into SBOL compliant format
    Each SBOLdocument is comprised of only one pathway
    """
    dict_reaction = None
    dict_chemicals = None

    def __init__(self, rxn_dat_file, chem_dat_file):
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

    def add_new_path(self, reactions_str, pathway_param):
        lines = reactions_str.split('\n')
        index = pathway_param[0]

        namespace = 'http://sbolstandard.org/testfiles' + str(index)
        sbol3.set_namespace(namespace)
        doc = sbol3.Document()
        intro = sbol3.TopLevel('pathway_properties', namespace + '/' + 'pathway_analysis')
        doc.add(intro)
        i = 0
        for line in lines:
            if line == '':
                continue
            self.add_rxn(line, doc, namespace, i)
            i += i

        doc.write('pathway' + str(index) + '.xml')

    def add_rxn(self, rxn_str, SBOL_doc, namespace, i):
        line = rxn_str.split('\t')
        reaction = line[0]
        reaction_str = reaction.replace('-', '_')
        reaction_component = sbol3.Component(reaction_str, sbol3.SBO_BIOCHEMICAL_REACTION)
        reaction_component.description = self.dict_reaction.get(reaction)
        SBOL_doc.add(reaction_component)
        interaction = sbol3.Interaction('rxn' + str(i))
        interaction.type = sbol3.SBO_BIOCHEMICAL_REACTION
        chemicals = line[1].split('-->')
        reactants = chemicals[0]
        products = chemicals[1]

        if type(reactants) == str:
            reactants = [reactants]

        if type(products) == str:
            products = [products]

        for reactant in reactants:
            reactant = reactant.replace(' ', '')
            reactant_str = reactant.replace('-', '_')
            if SBOL_doc.find(namespace + '/' + reactant_str) is not None:
                reactant_component = SBOL_doc.find(namespace + '/' + reactant_str)
            else:
                reactant_component = sbol3.Component(reactant_str, sbol3.SBO_SIMPLE_CHEMICAL)
                compound = self.dict_chemicals.get(reactant)
                reactant_component.sequences = [
                    sbol3.Sequence(reactant_str + 'seq', elements=compound.inchi, encoding=sbol3.INCHI_ENCODING)]
                SBOL_doc.add(reactant_component)
            interaction.participations = [sbol3.Participation(reactant_str + '_reactant', sbol3.SBO_REACTANT)]
            reactant_sc = sbol3.SubComponent(reactant_component)
            reaction_component.features += [reactant_sc]

        for product in products:
            product = product.replace(' ', '')
            product_str = product.replace('-', '_')
            if SBOL_doc.find(namespace + '/' + product_str) is not None:
                product_component = SBOL_doc.find(namespace + '/' + product_str)
            else:
                product_component = sbol3.Component(product_str, sbol3.SBO_SIMPLE_CHEMICAL)
                compound = self.dict_chemicals.get(product)
                product_component.sequences = [
                    sbol3.Sequence(product_str + 'seq', elements=compound.inchi, encoding=sbol3.INCHI_ENCODING)]

            interaction.participations = [sbol3.Participation(product_str + '_product', sbol3.SBO_PRODUCT)]
            product_sc = sbol3.SubComponent(product_component)
            reaction_component.features += [product_sc]

        reaction_component.interactions = [interaction]


if __name__ == "__main__":
    rxn_dat = '/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/reactions.txt'
    chem_dat = '/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/chems.txt'
    SBOLDocumenter = SBOLDocumenter(rxn_dat, chem_dat)
    SBOLDocumenter.add_new_path("RXN-161\tCPD-347 --> CPD-347\n", [1, 2, 3])
