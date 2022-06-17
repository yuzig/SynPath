class Chemicals:
    unique_id = ""
    common_name = ""
    inchi = ""
    smiles = ""

    def __init__(self, unique_id, common_name, inchi, smiles):
        self.unique_id = unique_id
        self.common_name = common_name
        self.inchi = inchi
        self.smiles = smiles
