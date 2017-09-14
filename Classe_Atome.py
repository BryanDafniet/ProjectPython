#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###########################
# importation des modules #
###########################

import sys
import os

#######################
#Définition des objets#
#######################

class Residus(object):
	def __init__(self,resid_name,resid_num, list_atom):
		'''Description '''
		self.residu_name = resid_name
		self.num = resid_num
		self.tab_atom = []
	 	# boucle de lecture
		for i in range(len(list_atom)):
			#condition
			self.tab_atom.append(Atoms(list_atom[i][0], list_atom[i][1],
			list_atom[i][2], list_atom[i][3]))

class Atoms():
	def __init__(self, atom_name ,coord_X ,coord_Y , coord_Z ):
		# Residus.__init__(self,resid_name,resid_num)
		# self.res_name = resid_name
		# self.res_num = resid_num
		self.X = coord_X
		self.Y = coord_Y
		self.Z = coord_Z
		self.atom_name = atom_name

class Lecture_pdb:
	def __init__(self,file_src):
		self.tab_res = []
		with open(file_src,'r') as src :
			res_id = 0
			i = -1
			for line in src :
				if line[0:6].strip() in ("ATOM","TER") :
					if (res_id != int(line[22:26].strip()) or
					 	line[0:6].strip() == "TER") :
						if 'resid_name' in locals() and len(list_atom) == 4:
							print (list_atom)
							self.tab_res.append(
										Residus(resid_name,resid_num,list_atom))
						list_atom = {}
						j = 0
						i += 1
						res_id = int(line[22:26].strip())
					if (line[12:16].strip() in ("C", "H", "N", "O") and
						line[16:17] in (" ","A")): #Prise en compte uniquement
												   #de la position alternative A
						list_atom[j] = ([line[12:16].strip(),
										float(line[30:38].strip()),
										float(line[38:46].strip()),
										float(line[46:54].strip())])
						j += 1
						resid_name = line[17:20]
						resid_num = res_id

	def extract_tab(self):

		return self.tab_res



class Check_and_prepare_pdb(object):
	""" Cette classe permet de valider la ligne de commande entrée.
	Ainsi que de confirmer la présence et ajouter les hydrogènes au fichier pdb.
	Enfin il initialise la lecture du fichier.
	"""
	def __init__(self,pdb_file):

		if len(sys.argv) != 2:
			print ("Usage: ./Classe_Atome.py fichier.pdb")
			sys.exit()
		self.file_name = pdb_file
		self.new_file = pdb_file[:-4] + "_H.pdb"

		if not os.path.exists(self.file_name):
			sys.exit("Le fichier " + self.file_name + " est absent")
		os.system("./reduce.3.23.130521 " + self.file_name + " >"
				  + self.new_file)#Ajout Hydrogènes

	def lecture(self):

		return Lecture_pdb(self.new_file).extract_tab()

######################
# Création fonctions #
######################

# def write_pdb_chain(liste_chain,file_dest,count):
# 	if count == 0 :
# 		with open(file_dest,'w') as dest:
# 			for i in xrange(len(liste_chain)):
# 				dest.write(liste_chain[i])
# 	elif count != 0:
# 		with open(file_dest,'a') as dest:
# 			for i in xrange(len(liste_chain)):
# 				dest.write(liste_chain[i])


class dist_energy(Residus):
	def __init__()
	def distAtom(self, resid_name, tabAtom):
		for i in range(0,len(Atoms)-1): #i = iterateur nombre residus
			for j in range(i+1,len(Atoms)): # j = iterateur nombre residus à comparer
				dist_OH = math.sqrt(((Atoms[j].tabAtom[0].X)- #O = 0 H = 1
				(Atoms[i].tabAtom[1].X))**2 + ((Atoms[j].tabAtom[0].Y)-
				(Atoms[i].tabAtom[1].Y))**2 + ((Atoms[j].tabAtom[0].Z)-
				(Atoms[i].tabAtom[1].Z))**2)
				if dist_OH <= 3.5 and (Atoms[i] and Atoms[j] != "PRO"):
					dist_ON = math.sqrt(((Atoms[j].tabAtom[0].X)- #O = 0 N = 2
					(Atoms[i].tabAtom[2].X))**2 + ((Atoms[j].tabAtom[0].Y)-
					(Atoms[i].tabAtom[2].Y))**2 + ((Atoms[j].tabAtom[0].Z)-
					(Atoms[i].tabAtom[2].Z))**2)

					dist_CH = math.sqrt(((Atoms[j].tabAtom[3].X)- #C = 3 H = 1
					(Atoms[i].tabAtom[1].X))**2 + ((Atoms[j].tabAtom[3].Y)-
					(Atoms[i].tabAtom[1].Y))**2 + ((Atoms[j].tabAtom[3].Z)-
					(Atoms[i].tabAtom[1].Z))**2)

					dist_CN = math.sqrt(((Atoms[j].tabAtom[3].X)- #C = 3 N = 2
					(Atoms[i].tabAtom[2].X))**2 + ((Atoms[j].tabAtom[3].Y)-
					(Atoms[i].tabAtom[2].Y))**2 +	((Atoms[j].tabAtom[3].Z)-
					(Atoms[i].tabAtom[2].Z))**2)
					return Atoms[i],Atoms[j],dist_OH, dist_CN,dist_CH,dist_ON
				else:
					break

	def	energie(self,Res_1,Res_2,dist_OH,dist_CN,dist_CH,dist_ON):
		energie = 0.084*((1/dist_ON)+(1/dist_CH)-(1/dist_OH)-1/dist_CN)*332
		if energie <= -0.5:
			print ("Il y'a un H-bond entre {} et {}".format(Res_1,Res_2))
		else:
			print ("Pas de H-bond entre {} et {}".format(Res_1, Res2))

##################
# prog principal #
##################

x=Check_and_prepare_pdb(sys.argv[1]).lecture()
for i in range(len(x)):
	print(x[i].residu_name)
	for j in range(4):
			 	print(x[i].tab_atom[j].atom_name,x[i].tab_atom[j].X,x[i].tab_atom[j].Y,x[i].tab_atom[j].Z)

# test=Lecture_pdb("2am9_H.pdb").extract_tab()
# print(test[1].residu_name)

# if len(sys.argv) != 2:
# 	print_usage_and_quit()

# check_file(sys.argv[1])
# file_name = sys.argv[1]
# chain_id=[]
# chain_pdb=[]
# chain_id_name=""
# count=0
# for i in xrange(2,len(sys.argv)):
# 	chain_id.append(sys.argv[i].upper())

# for j in chain_id:
# 	chain_id_name += "_"+j

# for j in chain_id:
# 	chain_pdb=extract_chain(file_name,j,file_name[:-4]+chain_id_name+".pdb")
# 	write_pdb_chain(chain_pdb,file_name[:-4]+chain_id_name+".pdb",count)
# 	count += 1
