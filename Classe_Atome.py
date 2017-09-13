#!/usr/bin/env python3

class Residus(object):
	def __init__(self,resid_name,resid_num,list_atom):
		'''Description '''
		self.residu_name = resid_name
		self.num = resid_num
		self.tabAtom = []
	# def lecture_pdb(list_atom):
	# 	# boucle de lecture
		for i in range(4):
			#condition
			self.tabAtom.append(Atoms(list_atom[i][0],list_atom[i][1],
			list_atom[i][2], list_atom[i][3]))


class Atoms(Residus):
	def __init__(self, atom_name ,coord_X ,coord_Y , coord_Z ):
		# Residus.__init__(self,resid_name,resid_num)
		# self.res_name = resid_name
		# self.res_num = resid_num
		self.X = coord_X
		self.Y = coord_Y
		self.Z = coord_Z
		self.atom_name = atom_name

atom={}
atom[0] = ["CA",1,2,3]
atom[1] = ["CA",1,2,3]
atom[2] = ["N",-2,0,1]
atom[3] = ["C",1,-3,2]
coord_Y=atom[1][2]
coord_Z=atom[1][3]
coord_X=atom[1][1]
resid_num=1
residu_1=Residus("ILE",2, atom)
print (residu_1.tabAtom[0].atom_name)
# test=Atoms('CA',1,2,3)
# test=test+Atoms("TRP",1,'N',-1,-3,2)
# print(test.residu_name)
# for i in test,test1:
# 	if i.res_name=="TRP":
# 		print (i.res_name,i.atom_name)
