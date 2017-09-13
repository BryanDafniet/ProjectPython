#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from math import sqrt

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


	def distAtom(self,Res):
		for i in range(0,len(Res)-1): #i = iterateur nombre residus
			for j in range(i+1,len(Res)): # j = iterateur nombre residus à comparer
				dist_OH = math.sqrt(((Res[j].tabAtom[0].X)- #O = 0 H = 1
				(Res[i].tabAtom[1].X))**2 + ((Res[j].tabAtom[0].Y)-
				(Res[i].tabAtom[1].Y))**2 + ((Res[j].tabAtom[0].Z)-
				(Res[i].tabAtom[1].Z))**2)
				if dist_OH <= 3.5 and (Res[i] and Res[j] != "PRO"):
					dist_ON = math.sqrt(((Res[j].tabAtom[0].X)- #O = 0 N = 2
					(Res[i].tabAtom[2].X))**2 + ((Res[j].tabAtom[0].Y)-
					(Res[i].tabAtom[2].Y))**2 + ((Res[j].tabAtom[0].Z)-
					(Res[i].tabAtom[2].Z))**2)	

					dist_CH = math.sqrt(((Res[j].tabAtom[3].X)- #C = 3 H = 1
					(Res[i].tabAtom[1].X))**2 + ((Res[j].tabAtom[3].Y)-
					(Res[i].tabAtom[1].Y))**2 + ((Res[j].tabAtom[3].Z)-
					(Res[i].tabAtom[1].Z))**2)

					dist_CN = math.sqrt(((Res[j].tabAtom[3].X)- #C = 3 N = 2
					(Res[i].tabAtom[2].X))**2 + ((Res[j].tabAtom[3].Y)-
					(Res[i].tabAtom[2].Y))**2 +	((Res[j].tabAtom[3].Z)-
					(Res[i].tabAtom[2].Z))**2)
					return Res[i],Res[j],dist_OH, dist_CN,dist_CH,dist_ON
				else:
					break

	def	energie(self,Res_1,Res_2,dist_OH,dist_CN,dist_CH,dist_ON):
		energie = 0.084*((1/dist_ON)+(1/dist_CH)-(1/dist_OH)-1/dist_CN)*332
		if energie < -0.5:
			print ("Il y'a un H-bond entre {} et {}".format(Res_1,Res_2))
		else:
			print ("Pas de H-bond entre {} et {}".format(Res_1, Res2))


atom1 = {}
atom2 = {}
atom1[1] = ["H",1,2,3]
atom1[0] = ["O",1,2,3]
atom1[2] = ["N",-2,0,1]
atom1[3] = ["C",1,-3,2]

residu_1=Residus("ILE",2, atom1)

# test=Atoms('CA',1,2,3)
# test=test+Atoms("TRP",1,'N',-1,-3,2)
# print(test.residu_name)
# for i in test,test1:
# 	if i.res_name=="TRP":
# 		print (i.res_name,i.atom_name)
atom2[1] = ["H",1.1,2.2,3.3]
atom2[0] = ["O",1.2,2.5,3.4]
atom2[2] = ["N",-2.3,0.1,1.2]
atom2[3] = ["C",1.4,-3.5,2.6]

print Res
residu_2=Residus("MET",3, atom2)
Res = [residu_1,residu_2]
print (residu_1.tabAtom[0].atom_name)
print (Res.distAtom)