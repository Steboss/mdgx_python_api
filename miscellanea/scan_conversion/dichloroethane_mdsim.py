#script for ing to remember:popen may be true when it's finished, ratgher than false.
#preparing input file of dihchloroethane/forcefield project. Go into the md folder with dichlorethane.mol2 and box creation folder already prepared  and  type python dicghloroethane_mdsim.py  The script use tleap for system creation, loading the new frcmod file. Then the box is preparaed with sander according to FESetup 
#one good th
import os, shutil, sys, re, glob, threading, time
import subprocess as subp



def sanderminmd(folder):
    box_folder=folder+"/box_creation"
    command = "sander -i min00001.in -p solvated.parm7 -c solvated.rst7 -O -o min00001.out -e min00001.en -x min00001.nc -inf min00001.info -r min00001.rst7 -ref solvated.rst7"
    print("Running sander minimization step 1 out of 4")
    proc= subp.Popen(command, cwd=box_folder, shell=True)
    proc.wait()
     
###########heat the system up############
    command = "sander -i md00002.in -p solvated.parm7 -c min00001.rst7 -O -o md00002.out -e md00002.en -x md00002.nc -inf md00002.info -r md00002.rst7 -ref solvated.rst7"
    print("Running sander heat system step 2 out of 4")
    proc = subp.Popen(command, cwd=box_folder, shell=True)
    proc.wait()
   
###########md############ 
    command = "sander -i md00003.in -p solvated.parm7 -c md00002.rst7 -O -o md00003.out -e md00003.en -x md00003.nc -inf md00003.info -r md00003.rst7 -ref solvated.rst7"
    print("Running sander heat system step 3 out of 4")
    proc = subp.Popen(command, cwd=box_folder, shell=True)
    proc.wait()
   
###########pressurize############
    command = "sander -i md00004.in -p solvated.parm7 -c md00003.rst7 -O -o md00004.out -e md00004.en -x md00004.nc -inf md00004.info -r md00004.rst7 -ref md00003.rst7"
    print("Running sander md step 4 out of 4")
    proc = subp.Popen(command, cwd=box_folder, shell=True)
    proc.wait()
    
if __name__ == "__main__":
 
    mol2_file = "dichloroethane.mol2"
    tleap_dat = open("tleap.dat", "w")
    
    files = os.listdir( os.getcwd() )
    files.sort()

    for f in files:
        if f.endswith(".frcmod"):
            name, ftype = f.split(".frcmod")
            print("Found frcmod file %s \n" % name)
            frcmod_file = name + ".frcmod"
        
    print("tleap: creation of box for file %s " % mol2_file)
    
    tleap_dat.write(" source leaprc.gaff\n source leaprc.ff14SB\n a = loadmol2 %s\n laodamberparams %s\n saveamberparm a dichloroethane.prmtop dichloroethane.inpcrd\n loadamberparams frcmod.ionsjc_tip3p\n loadamberparams frcmod.ionslrcm_hfe_tip3p\n solvatebox a TIP3PBOX 12.0 0.75 \n addIons a Na+ 0 \n addIons a Cl- 0 \n saveamberparm a solvated.parm7 solvated.rst7\n quit\n" %(mol2_file,frcmod_file))
           
    command = "tleap -f tleap.dat"
    print("%s" % command)
    current_folder = os.getcwd()
    proc = subp.Popen(command,cwd=current_folder, shell=True, stdout = subp.PIPE, stderr = subp.PIPE)
    if  proc:
        print("Copying file for sander")
        shutil.move("solvated.parm7","box_creation")
        shutil.move("solvated.rst7","box_creation")
        sanderminmd(current_folder)
      #tleapcreator(current_folder,command)
      #sanderminmd(current_folder)
   
     
