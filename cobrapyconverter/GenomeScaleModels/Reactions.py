class Reactions:
    ecnum = ''
    reactants = ''
    products = ''
    unique_id = ''
    direction = ''

    def __init__(self, ecnum, reactants, products, unique_id, direction):
        self.ecnum = ecnum
        self.reactants = reactants
        self.products = products
        self.unique_id = unique_id
        self.direction = direction
