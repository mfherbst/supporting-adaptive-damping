using PyCall
super = pyimport("ase.build").bulk("Cr", cubic=true) * (1, 1, 10)
super.pop(10)
super.write("Cr19.json")
