#!/usr/bin/env python3

class Residus(object):
	'''Description du résidus'''
	def __init__(self, nameAA = "MET", n = 450):
		Residus.nameAA = nameAA
		self.number = n
	def show():
		print ("Nom Residus:", Residus.nameAA, "\n")
	

	class Atom():
		'''Définition des différentes caractéristiques de l'atome'''
		def __init__(self, atom_name = "CA", coord_X = 0.1, coord_Y = 0, coord_Z = 0.3):
			Residus.__init__(self)
			self.name = atom_name
			self.X = coord_X
			self.Y = coord_Y
			self.Z = coord_Z
		def show(self):
			Residus.show()
			print ("Type d'atome :", self.name, "\nCoordonnée X :",
			 self.X, "\nCoordonnée Y :", self.Y, "\nCoordonnée Z :", self.Z)
			
# test = Atom()
# print (test.nameAA)
c1 = Residus.Atom(atom_name = "TEST")
c1.show()

class distAtom():
	Residus
