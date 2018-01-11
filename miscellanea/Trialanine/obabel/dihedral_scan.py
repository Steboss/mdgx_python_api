import os,re, sys, shutil
import math
import numpy

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *

from Sire.Tools.DCDFile import *

from Sire.Tools import Parameter, resolveParameters

import Sire.Stream



def createSystem(molecules):

    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system


def generateDihedrals(old_solute,atom_list,atom_list_2): #residue_list,atom_list):

    connectivity = old_solute.property("connectivity")
    all_dihedrals = connectivity.getDihedrals()
    atom0 = atom_list[0]
    atom1 = atom_list[1]
    atom2 = atom_list[2]
    atom3 = atom_list[3]

    atom4 = atom_list_2[0]
    atom5 = atom_list_2[1]
    atom6 = atom_list_2[2]
    atom7 = atom_list_2[3]




    for dihedral in all_dihedrals:
#        print(dihedral)
        if (atom4 == dihedral.atom0() and atom5 == dihedral.atom1() and atom6 == dihedral.atom2() and atom7 == dihedral.atom3()):
            print("OK")
            dihbond = BondID(dihedral.atom1(), dihedral.atom2())
            for deg in range(180,540,10): #here we should put the initial value of the dihedral
                new_solute = old_solute.move().change(dihbond,deg*degrees).commit()
                #Here for each new_solute  you have to restrain and minimized and then save:
                #Apply restraint   1) SpecRest
                #Minimization      2) Minimization
                #Save              3) GenerateCoordFiles
                generateCoordFiles(new_solute.property("coordinates"),deg)


def generateCoordFiles(coordinates,idx):


     coord_toVect = coordinates.toVector()
     string_coord = str(coord_toVect)
     string_rev   = re.sub("[^0-9.e-]", " ", string_coord)
     string_split = string_rev.split()


     single_idx = idx
     file_name = "amber_inps/alanine" + str(single_idx) + ".mdcrd"
     singlecoord = open(file_name, "w")
     singlecoord.write("ACE\n")
     singlecoord.write("    42\n")
     index=0
     for f in string_split:
         if float(f)<0:
             singlecoord.write("  %.7f" % float(f))
         elif float(f)>=10:
             singlecoord.write("  %.7f" % float(f))
         else:
             singlecoord.write("   %.7f" % float(f))
         index+=1
         if not index%6:
             singlecoord.write("\n")



##########
###MAIN###
##########

if __name__ == "__main__":


     try:
         top_file = sys.argv[1]
         crd_file = sys.argv[2]


     except IndexError:

         top_file = "SYSTEM.top"
         crd_file = "SYSTEM.crd"


     amber = Amber()

     molecules, space = amber.readCrdTop(crd_file, top_file)
     system = createSystem(molecules)
     old_solute = system[MGName("all")].moleculeAt(0).molecule()
     residues = old_solute.residues()
     print("This are the atoms")
     for res in residues:

         print(res)
         for atom in res.atoms():
             print(atom)

     accept=True
     print("Now give the residue-atom in this way:")
     print("Example: ALA2-CA,ALA3-N,ALA3-CA,ALA3-C")
     while(accept):

         phi = input("What are the residue-atom for PHI?")
         psi = input("What are the residue-atom for PSI")
         print("Here the selection:")
         print("PHI:")
         print(phi)
         print("Psi:")
         print(psi)
         confirm = input("Is it correct? yes/no\n")
         if confirm=="yes":
             accept=False
         else:
             accept=True #so we can refix the angle


     #what we need is
     #1) take as input phi and psi
     #2) understand which angles move more atoms and keep it freeze
     #3) make the scan of the "less" atoms angles
     #4) obtain phi and psi = 0 or take the initial value and scan them form there

     atom_list = []

     if residues.count() != 1:  #means that there are residues
        for dihedral_atom in (dihedral_list.split(",")):#Gave as ALA2-CA,ALA3-N,ALA3-CA,ALA3-CB
             residue = str(dihedral_atom.split("-")[0])
             #print(residue)
             res_numb = str((re.findall(r'(\d+)',residue))[0])
             #print(res_numb) #this is the residue number
             res_name = residue.split(res_numb)[0]
             res_sire = res_name + " : " + res_numb
             #print(res_sire)
             atom = dihedral_atom.split("-")[1]
             for i in range(0, residues.count()):
                 res_value = str(residues[i].name().value()) + " : " + str(residues[i].number().value())  # residue name from Sire
                 if res_sire == res_value:
                     atoms = residues[i].atoms()
                     for at in range(0, atoms.count()):
                         if atom == atoms[at].name().value():
                             atom_idx  = residues[i].atomIdxs()[at]
                             atom_list.append(atom_idx)


     print(atom_list)
     #generateDihedrals(old_solute,atom_list,bigfile)

     choice=input("Would like to scan another dihedral?")
     if choice=="yes":
         dihedral_list = input("What dihedral do you want to scan?")    #let's begin with just one dihedral for a test be carefule raw_input = input  in new Python
         atom_list_2 = []
         if residues.count() != 1:  #means that there are residues
            for dihedral_atom in (dihedral_list.split(",")):#Gave as ALA2-CA,ALA3-N,ALA3-CA,ALA3-CB
                residue = str(dihedral_atom.split("-")[0])
                #print(residue)
                res_numb = str((re.findall(r'(\d+)',residue))[0])
                #print(res_numb) #this is the residue number
                res_name = residue.split(res_numb)[0]
                res_sire = res_name + " : " + res_numb
                #print(res_sire)
                atom = dihedral_atom.split("-")[1]
                for i in range(0, residues.count()):
                    res_value = str(residues[i].name().value()) + " : " + str(residues[i].number().value())  # residue name from Sire
                    if res_sire == res_value:
                        atoms = residues[i].atoms()
                        for at in range(0, atoms.count()):
                            if atom == atoms[at].name().value():
                                atom_idx  = residues[i].atomIdxs()[at]
                                atom_list_2.append(atom_idx)


     print(atom_list_2)

     generateDihedrals(old_solute,atom_list,atom_list_2)
