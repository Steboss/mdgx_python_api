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


def generateDihedrals(old_solute,atom_list,bigfile): #residue_list,atom_list):

    connectivity = old_solute.property("connectivity")
    all_dihedrals = connectivity.getDihedrals()
    atom0 = atom_list[0]
    atom1 = atom_list[1]
    atom2 = atom_list[2]
    atom3 = atom_list[3]


    
    
    for dihedral in all_dihedrals:
        if (atom0 == dihedral.atom0() and atom1 == dihedral.atom1() and atom2 == dihedral.atom2() and atom3 == dihedral.atom3()):
            print("OK")
            dihbond = BondID(dihedral.atom1(), dihedral.atom2())
            for deg in range(0,360):
                new_solute = old_solute.move().change(dihbond,deg*degrees).commit()
                #Here for each new_solute  you have to restrain and minimized and then save:
                #Apply restraint   1) SpecRest
                #Minimization      2) Minimization
                #Save              3) GenerateCoordFiles
                generateCoordFiles(bigfile,new_solute.property("coordinates"),deg)
                
                
def generateCoordFiles(bigfile,coordinates,idx):
     
     
     index = 0
     
     coord_toVect = coordinates.toVector()
     string_coord = str(coord_toVect)
     string_rev   = re.sub("[^0-9.e-]", " ", string_coord)
     string_split = string_rev.split()
     
     for f in string_split:
         if float(f)<0:
             bigfile.write("  %.3f" % float(f))
         elif float(f)>=10:
             bigfile.write("  %.3f" % float(f))
         else:
             bigfile.write("   %.3f" % float(f))
         index+=1
         if not index%10:
             bigfile.write("\n")
     bigfile.write("\n")
        
     single_idx = idx
     file_name = "amber_inps/alanine" + str(single_idx) + ".mdcrd"
     singlecoord = open(file_name, "w")
     singlecoord.write("LIG\n")
     singlecoord.write("    32\n")
     for f in string_split:
         if float(f)<0:
             singlecoord.write("  %.6f" % float(f))
         elif float(f)>=10:
             singlecoord.write("  %.6f" % float(f))
         else:
             singlecoord.write("   %.6f" % float(f))
         index+=1
         if not index%6:
             singlecoord.write("\n")
 
########################################MAIN#################################################################

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

     dihedral_list = input("What dihedral do you want to scan?")    #let's begin with just one dihedral for a test be carefule raw_input = input  in new Python
    

     residues = old_solute.residues()
     
     atom_list = []
     bigfile = open("ALANINE.mdcrd", "a")
     bigfile.write("Cpptraj generator\n")
  
     if residues.count() != 1:  #means that there are residues
         for dihedral_atom in (dihedral_list.split()):
             for i in range(0, residues.count()):
                 splitter = dihedral_atom.split(":")
                 residue = str(splitter[0])
                 atom = splitter[1]
                 res_numb_split = re.split("(\d+)", residue)  #division of residue name and residue number
                 res_name = res_numb_split[0] + " : " + res_numb_split[1] # creation of Sire residue name e.g. ( ALA : 2 ) 
                 res_value = str(residues[i].name().value()) + " : " + str(residues[i].number().value())  # residue name from Sire
                 
                 
                 if res_name == res_value:
                     atoms = residues[i].atoms()
                     for at in range(0, atoms.count()):
                         if atom == atoms[at].name().value():
                             atom_idx  = residues[i].atomIdxs()[at]
                             atom_list.append(atom_idx)
     
             
     generateDihedrals(old_solute,atom_list,bigfile)   





        #getdihedral: system, residue list, atom list

#This will be useful        if use_restraints.val:
                                 #system = setupRestraints(system)
        
      #  if use_spec_restraint.val:
       #     system = specificRestraint(system)
        #    print("System has been restrained")
     
