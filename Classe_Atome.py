#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###########################
# importation des modules #
###########################

import sys
import os
import math

#########################
# Définition des objets #
#########################


class Residus(object):
    def __init__(self, resid_name, resid_num, list_atom):
        '''Description '''
        self.residu_name = resid_name
        self.num = resid_num
        self.tab_atom = []
        for i in range(len(list_atom)):
            self.tab_atom.append(Atoms(list_atom[i][0], list_atom[i][1],
                                       list_atom[i][2], list_atom[i][3]))


class Atoms():
    def __init__(self, atom_name, coord_X, coord_Y, coord_Z):
        self.X = coord_X
        self.Y = coord_Y
        self.Z = coord_Z
        self.atom_name = atom_name

    def dist_atome(self, atom_2):
        dist = math.sqrt((self.X - atom_2.X)**2 +
                         (self.Y - atom_2.Y) ** 2 +
                         (self.Z - atom_2.Z) ** 2)
        return dist


class Lecture_pdb:
    def __init__(self, file_src):
        self.tab_res = []
        with open(file_src, 'r') as src:
            res_id = 0
            i = -1
            for line in src:
                if line[0:6].strip() in ("ATOM", "TER"):
                    if(res_id != int(line[22:26].strip()) or
                       line[0:6].strip() == "TER"):
                        if 'resid_name' in locals() and len(list_atom) == 4:
                            self.tab_res.append(
                                        Residus(resid_name,
                                                resid_num,
                                                list_atom))
                        list_atom = {}
                        j = 0
                        i += 1
                        res_id = int(line[22:26].strip())
                    # Prise en compte uniquement de la position alternative A
                    if(line[12:16].strip() in ("C", "H", "N", "O") and
                       line[16:17] in (" ", "A")):
                        list_atom[j] = ([line[12:16].strip(),
                                        float(line[30:38].strip()),
                                        float(line[38:46].strip()),
                                        float(line[46:54].strip())])
                        j += 1

                        resid_name = translate_code(line[17:20])
                        resid_num = res_id

    def extract_tab(self):

        return self.tab_res


class Check_and_prepare_pdb(object):
    """ Cette classe permet de valider la ligne de commande entrée.
    Ainsi que de confirmer la présence et ajouter les hydrogènes au fichier
    pdb. Enfin il initialise la lecture du fichier.
    """
    def __init__(self, pdb_file):

        if len(sys.argv) != 2:
            print("Usage:./Classe_Atome.py fichier.pdb")
            sys.exit()
        self.file_name = pdb_file
        self.new_file = pdb_file[:-4] + "_H.pdb"

        if not os.path.exists(self.file_name):
            sys.exit("Le fichier " + self.file_name + " est absent")
        os.system("./reduce.3.23.130521 " + self.file_name +
                  " >" + self.new_file)  # Ajout Hydrogènes

    def lecture(self):

        return Lecture_pdb(self.new_file).extract_tab()


######################
# Création fonctions #
######################
def translate_code(code_3_lettre):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    code_1_lettre = d[code_3_lettre]
    return code_1_lettre


def calcule_distance(list_res):
    list_distance = []
    # i = iterateur nombre residus
    # j = iterateur nombre residus à comparer
    for i in range(len(list_res) - 1):
        for j in range(i + 1, len(list_res)):
            dist_OH = (list_res[j].tab_atom[3].
                       dist_atome(list_res[i].tab_atom[2]))
            if dist_OH <= 3.5:
                dist_ON = (list_res[j].tab_atom[0].
                           dist_atome(list_res[i].tab_atom[2]))

                dist_CH = (list_res[j].tab_atom[3].
                           dist_atome(list_res[i].tab_atom[1]))

                dist_CN = (list_res[j].tab_atom[0].
                           dist_atome(list_res[i].tab_atom[1]))
                dico_dist = {}
                dico_dist['res_num1'] = list_res[i].num
                dico_dist['res_num2'] = list_res[j].num
                dico_dist['dist_OH'] = dist_OH
                dico_dist['dist_CN'] = dist_CN
                dico_dist['dist_CH'] = dist_CH
                dico_dist['dist_ON'] = dist_ON
                dico_dist['res_name1'] = list_res[i].residu_name
                dico_dist['res_name2'] = list_res[j].residu_name
                list_distance.append(dico_dist)
    return list_distance


def energie(list_distance):
    list_hbond = []
    for distance in list_distance:
        energie = 0.084 * ((1 / distance["dist_ON"]) +
                           (1 / distance["dist_CH"]) -
                           (1 / distance["dist_OH"]) -
                           (1 / distance["dist_CN"])) * 332
        if energie <= -0.5:
            list_hbond.append(distance)
    return list_hbond


def helix(list_hbond):
    list_helix=[]
    for i in range(len(list_hbond)):  # parcours de la liste
        for j in (3,4,5):  # helix 3, 4 ou 5
            if list_hbond[i]["res_num1"] == list_hbond[i]["res_num2"] - j:
                k= i+1
                print("test",i)
                while (list_hbond[k]["res_num1"] <
                       list_hbond[i]["res_num1"] + 2):
                    if ((list_hbond[k]["res_num1"] ==
                        list_hbond[i]["res_num1"] + 1) and
                        (list_hbond[k]["res_num1"] ==
                         list_hbond[k]["res_num2"] - j)):
                        for l in (i,k):
                            dico_helix = {}
                            dico_helix['res_num'] = list_hbond[l]["res_num1"]
                            dico_helix['res_name'] = list_hbond[l]["res_name1"]
                            dico_helix["turn"] = j
                            list_helix.append(dico_helix)

                    k += 1
                    i = k

                print("ta_maman",i)




##################
# prog principal #
##################

list_atoms = Check_and_prepare_pdb(sys.argv[1]).lecture()
list_dist = calcule_distance(list_atoms)
list_hbond = energie(list_dist)
print(len(list_dist))
print(list_hbond[1])
helix(list_hbond)
# for i in range(len(x)):
#     print(x[i].residu_name, x[i].num)
#     for j in range(4):
#                 print(x[i].tab_atom[j].atom_name, x[i].tab_atom[j].X, x[i].tab_atom[j].Y, x[i].tab_atom[j].Z)
