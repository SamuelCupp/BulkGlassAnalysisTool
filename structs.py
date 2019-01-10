class atom:
    "Class for storing individual atoms"

    def __init__(self,element,location,VanderWaals):
        self.elem = element
        self.loc = location
        self.vdw = VanderWaals

class fragment:
    "Class for storing individual atoms"

    def __init__(self,fragName,offset):
        self.name = fragName
        self.offset = offset
        self.atomList = []

    def addAtom(self,element,location,VanderWaals):
        self.atomList.append(atom(element,location,VanderWaals))
