#!/usr/bin/env python3

class Residus(object):
	def __init__(self,resid_name = "ASP",resid_num = 0):
		'''Description '''
		self.resid_name = resid_name
		self.num = resid_num

class Atoms(Residus):
	def __init__(self,resid_name = "ASP",resid_num = 0,atom_name = "CA" ,coord_X = 0.1 ,coord_Y = 0.2, coord_Z = 0.3):
 		Residus.__init__(self,resid_name,resid_num)
 		self.X = coord_X
 		self.Y = coord_Y
 		self.Z = coord_Z
 		self.atom_name = atom_name

test=Atoms()
print (test.resid_name,test.atom_name)
