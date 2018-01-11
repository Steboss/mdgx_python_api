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

####Remember "to doctstring" everything, please#####

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



def atomMatching(oldsolute):

     atom_dictionary = {}  # to be functionalised
     atoms = old_solute.atoms()
     for i in range(0, atoms.count()):
         at = atoms[i]
         at_idx = at.index()
         at_name = at.name().value()
         at_number = int(at.number().value()) -1  #-1 due to mismatcing between chemistry.amber and sire
         atom_dictionary[str(at_name)] = [at_number] 
     return atom_dictionary


def atomAllMatching(oldsolute):

     atom_dictionary = {}  # to be functionalised
     atoms = old_solute.atoms()
     for i in range(0, atoms.count()):
         at = atoms[i]
         at_idx = at.index()
         at_name = at.name().value()
         at_number = int(at.number().value()) -1  #-1 due to mismatcing between chemistry.amber and sire
         atom_dictionary[str(at_idx)] = [ at_name, at_number]
     return atom_dictionary

def atomResidueCheck(oldsolute,dihedrallist):
     
     residues = oldsolute.residues()   #check over residues presence
     atom_list = []
     if dihedrallist == "all" or dihedrallist=="":   #check over dihedral presence
         atom_list = []
     else:
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
         else:
              for dihedral_atom in (dihedral_list.split()):
                  atom_list.append(dihedral_atom)
             
     return atom_list

def atomNameNumberList(dihedral_atom_list,atom_dictionary):


    atom_list_file = open("atom_list.dat","w")
    
    index=0
    
    for dihedral in dihedral_atom_list:
        vox = atom_dictionary[str(dihedral)]
        
        for v in vox:
            atom_list_file.write("%s    %s   " %(dihedral,str(v)))
            index+=1
            if index==8:
                atom_list_file.write("\n")
                index=0


def atomAllNameNumberList(dihedral_atom_list,atom_dictionary):

    atom_list_file = open("atom_list.dat","w")
    
    index=0
    
    for dihedral in dihedral_atom_list:
        vox = atom_dictionary[str(dihedral)]
        
        for v in vox:
            atom_list_file.write("%s   " %str(v))
            index+=1
            if index==8:
                atom_list_file.write("\n")
                index=0  



def dihedralcalculator(dihedral,old_solute):

    atm0 = dihedral.atom0()
    atm1 = dihedral.atom1()
    atm2 = dihedral.atom2()
    atm3 = dihedral.atom3()
    atm0 = old_solute.select(atm0)
    atm1 = old_solute.select(atm1)
    atm2 = old_solute.select(atm2)
    atm3 = old_solute.select(atm3)
    atm0coord = atm0.property("coordinates")
    atm1coord = atm1.property("coordinates")
    atm2coord = atm2.property("coordinates")
    atm3coord = atm3.property("coordinates")
    phi01 = space.calcDihedral(atm0coord,atm1coord,atm2coord,atm3coord)
    
    return math.ceil(phi01.to(degree))


def generateDihedrals(old_solute,bigfile,atom_dihedral,atom_dictionary): 
    connectivity = old_solute.property("connectivity")
    all_dihedrals = connectivity.getDihedrals()
    dihedral_atom_list = []

    if not atom_dihedral: #if the atom_dihedral list is empty perform a move on all dihedrals
        dihedral_len = len(all_dihedrals)
        index_order = 0
        for dihedral in all_dihedrals:
            atm0= dihedral.atom0()
            atm1= dihedral.atom1()
            atm2= dihedral.atom2()
            atm3= dihedral.atom3()
            dihedral_atom_list.append(atm0)
            dihedral_atom_list.append(atm1)
            dihedral_atom_list.append(atm2)
            dihedral_atom_list.append(atm3)
            dihbond = BondID(dihedral.atom1(), dihedral.atom2())
            for deg in range(0,180,10):     #scan every 10 degreees
                new_solute = old_solute.move().change(dihbond,deg*degrees).commit()
                generateCoordFiles(bigfile,new_solute.property("coordinates"),deg,index_order)
            index_order+=1
        atomAllNameNumberList(dihedral_atom_list,atom_dictionary)   
    else:
        index_dihedral = 0 
        index_order = 0 
        for numbdihedral in range(0, int(len(atom_dihedral)/4)):
            atm0 = atom_dihedral[index_dihedral]
            atm1 = atom_dihedral[index_dihedral+1]
            atm2 = atom_dihedral[index_dihedral+2]
            atm3 = atom_dihedral[index_dihedral+3]                               
            index_dihedral=index_dihedral+4
            
            dihedral_atom_list.append(atm0)
            dihedral_atom_list.append(atm1)
            dihedral_atom_list.append(atm2)
            dihedral_atom_list.append(atm3)
            for dihedral in all_dihedrals:
#                Find actual dihedral
#                print(old_solute.atom(dihedral.atom0()).name().value())
#                print(old_solute.atom(dihedral.atom1()).name().value())
#                print(old_solute.atom(dihedral.atom2()).name().value())
#                print(old_solute.atom(dihedral.atom3()).name().value())             
                if (atm0 == old_solute.atom(dihedral.atom0()).name().value() and atm1 == old_solute.atom(dihedral.atom1()).name().value() \
                   and atm2 == old_solute.atom(dihedral.atom2()).name().value() and atm3 == old_solute.atom(dihedral.atom3()).name().value()):
                    #Calculate the initial diehdral angle TODO: old_solute.select(at0) give Atom( Cl1:3 ) (example) could be useful and quick
                    initial_angle = dihedralcalculator(dihedral, old_solute)
                    dihbond = BondID(dihedral.atom1(), dihedral.atom2())
                    #print("initial dihedral %f " % initial_angle)
                    #print("final dihedral %f " % final_angle)
                    print("Found dihedral")
                    for deg in range(0,360,1):#initial_angle, final_angle,1):
                        actual_dihedral = initial_angle + deg
                        new_solute = old_solute.move().change(dihbond,deg*degrees).commit()  # you're adding degrees to your actual dihedral : (actual)initial_dihedral + deg
                        
                        generateCoordFiles(bigfile,new_solute.property("coordinates"),actual_dihedral,index_order)
                    index_order+=1
            atomNameNumberList(dihedral_atom_list,atom_dictionary) 
                
         
          
def generateCoordFiles(bigfile,coordinates,degrees,idx):
     
     natoms = coordinates.nAtoms()
     coord_toVect = coordinates.toVector()
     string_coord = str(coord_toVect)
     string_rev   = re.sub("[^0-9.e-]", " ", string_coord)
     string_split = string_rev.split()
     

     big_file_index = 0
     for f in string_split:
         if float(f)<0:
             bigfile.write("  %.3f" % float(f))
         elif float(f)>=10:
             bigfile.write("  %.3f" % float(f))
         else:
             bigfile.write("   %.3f" % float(f))
         big_file_index+=1
         if not big_file_index%10:
             bigfile.write("\n")
     bigfile.write("\n")

     small_file_index=0  
     
     if degrees >= 360:
         degrees = degrees - 360
      
     file_name = "amber_inps/ligand" + "_" + str(idx) + "_" + str(degrees) + ".mdcrd"   #for the moment let's call ligand the outputfile
    
     singlecoord = open(file_name, "w")
     singlecoord.write("LIG\n")
     singlecoord.write("    %s\n" % natoms) #coordinates
     for f in string_split:
         if float(f)<0:
             singlecoord.write("  %.7f" % float(f))   #absolutely important to have 7 decimals otherwise openmm minimizer does not work!!
         elif float(f)>=10:
             singlecoord.write("  %.7f" % float(f))
         else:
             singlecoord.write("   %.7f" % float(f))
         small_file_index+=1
         if not small_file_index%6:
             singlecoord.write("\n")

     #pdb_name = "amber_pdb_in/ligand" + "_" + str(idx) + "_" + str(degrees) + ".pdb"
     
     #pdbcoor = open(pdb_name, "w")

     #for f in string_split:
     #    pdbcoor.write("ATOM      1  C1  LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      2  C2  LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      3  Cl3 LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      4  Cl4 LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      5  H5  LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      6  H6  LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      7  H7  LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("ATOM      8  H8  LIG     1      %.3f   %.3f   %.3f  1.00  0.00\n" % float(f))
     #    pdbcoor.write("TER\n")
     #    pdbcoor.write("END\n")
  

 
########################################MAIN#################################################################

if __name__ == "__main__":
 

     try:
         top_file = sys.argv[1]     #it's better to give these input files if they are named  differently from SYSTEM.top/crd
         crd_file = sys.argv[2]
        
        
     except IndexError:
        
         top_file = "SYSTEM.top"
         crd_file = "SYSTEM.crd"
        

     output_crd_folder = "amber_inps"
     if not os.path.exists(output_crd_folder):
         os.makedirs(output_crd_folder)
     


     amber = Amber()
     
     molecules, space = amber.readCrdTop(crd_file, top_file)
     system = createSystem(molecules)
     old_solute = system[MGName("all")].moleculeAt(0).molecule()

     dihedral_list = input("What dihedral do you want to scan?")  # at the moment this option is not automatic. Future: automatic passing of the most influential dihedrals
     print("Modification of %s" % dihedral_list)

     bigfilename = str(top_file.split(".")[0]) + ".mdcrd"
     bigfile = open(bigfilename, "a")
     bigfile.write("Cpptraj generator\n")
     
     if dihedral_list == "all" or dihedral_list=="": 
         atom_dictionary = atomAllMatching(old_solute) #function for matching names of atoms between Sire and chemistry.amber
    
     else:
         atom_dictionary = atomMatching(old_solute)


     atom_dihedral = atomResidueCheck(old_solute,dihedral_list) #function for create atom list in presence of residues


     generateDihedrals(old_solute,bigfile,atom_dihedral,atom_dictionary) 

 
