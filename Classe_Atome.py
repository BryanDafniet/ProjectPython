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
				if line[0:6].strip() == "ATOM" :
					if res_id != int(line[22:26].strip()):
						if 'resid_name' in locals() and len(list_atom) == 4:
							self.tab_res.append(Residus(resid_name,resid_num,list_atom))
						list_atom = {}
						j = 0
						i += 1
						res_id = int(line[22:26].strip())
						if (line[12:16].strip() in ("C", "H", "N", "O") and
							line[16:17] in (" ","A")): #Prise en compte uniquement de la position alternative A
							list_atom[j] = ([line[12:16].strip(),
										float(line[30:38].strip()),
										float(line[38:46].strip()),
										float(line[46:54].strip())])
							j += 1
							resid_name =line[17:20]
							resid_num = res_id
							print (list_atom)
		# for i in range((len(self.tab_res))):
		# 	print(self.tab_res[i].tab_atom[0].atom_name)


class Check_and_prepare_pdb(object):
	""" Cette classe permet de valider la ligne de commande entrée.
	Ainsi que de confirmer la présence et ajouter les hydrogènes au fichier pdb.
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

		Lecture_pdb(self.new_file)








######################
# Création fonctions #
######################

# def extract_chain(file_src,id_chain,file_dest):
# 	with open(file_src,'r') as src:
# 		chain=[]
# 		for line in src:
# 			if	line[21:22].strip() == id_chain and ( line[0:6].strip() ==  "ATOM" or line[0:6].strip() == "TER" or line[0:6].strip() == "HETATM"):
# 				chain.append(line)
# 	if chain==[] :
# 		#os.remove(file_dest)
# 		sys.exit("La chaine "+id_chain+" n'est pas dans le fichier "+file_src)
# 	return chain

# def write_pdb_chain(liste_chain,file_dest,count):
# 	if count == 0 :
# 		with open(file_dest,'w') as dest:
# 			for i in xrange(len(liste_chain)):
# 				dest.write(liste_chain[i])
# 	elif count != 0:
# 		with open(file_dest,'a') as dest:
# 			for i in xrange(len(liste_chain)):
# 				dest.write(liste_chain[i])



##################
# prog principal #
##################

x=Check_and_prepare_pdb(sys.argv[1])
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
