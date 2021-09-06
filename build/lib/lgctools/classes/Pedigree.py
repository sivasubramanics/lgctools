
class Pedigree():

    def __init__(self, fone):
        self.name = fone

    def get_parent_a(self):
        return self.parent_a

    def get_parent_b(self):
        return self.parent_b

    def get_designation(self):
        return self.designation

    def set_parent_a(self, parent_a):
        self.parent_a = parent_a

    def set_parent_b(self, parent_b):
        self.parent_b = parent_b

    def set_designation(self, designation):
        self.designation = designation

    def __str__(self):
        return self.name
