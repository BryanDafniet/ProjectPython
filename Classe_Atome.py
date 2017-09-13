#!/usr/bin/env python3

class Residus(object):
	'''Description du résidus'''
	def __init__(self, nameAA = "MET", n = 450):
		Residus.nameAA = nameAA
		self.number = n
	def show():
		print ("Nom Residus:", Residus.nameAA, "\n")
	show = staticmethod(show)

class Atom(Residus):
	'''Définition des différentes caractéristiques de l'atome'''
	def __init__(self, atom_name = "CA", coord_X = 0.1, coord_Y = 0, coord_Z = 0.3):
		Residus.__init__(self)
		Atom.name = atom_name
		Atom.X = coord_X
		Atom.Y = coord_Y
		Atom.Z = coord_Z
	def show():
		Residus.show()
		print ("Type d'atome :", Atom.name, "\nCoordonnée X :",
		 Atom.X, "\nCoordonnée Y :", Atom.Y, "\nCoordonnée Z :", Atom.Z)
	show = staticmethod(show)	
# test = Atom()
# print (test.nameAA)
c1 = Atom()
c1.show()
