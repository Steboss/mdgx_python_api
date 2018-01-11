import sys,os,glob
import numpy



def write_submit(script,name):
    #name: the name of the structrue e.g. structure_180
    script.write("#!/bin/bash\n")
    script.write("#SBATCH -o output-%A-%a.out\n")
    script.write("#SBATCH -p serial -n 8\n")
    script.write("#SBATCH --time 48:00:00\n")
    script.write("\n")
    script.write("source /etc/profile.d/module.sh\n")
    script.write("export OMP_NUM_THREADS=1\n")
    script.write("g09 < %s.gzmat > %s.gout\n" %(name,name))
    script.write("\nwait\n")

def config_zmat(zmat,outputf,psi):

    outputf.write("%chk=structure.chk\n")
    outputf.write("NProcShared=8\n")
    outputf.write("#HF/6-31G* SCF=tight opt=Z-Matrix\n")
    outputf.write("\n")
    outputf.write("Scanning Phi 180:540 and Psi fixed at %d\n" %psi)
    outputf.write("\n")
    for lines in zmat[6:49]:
        outputf.write(lines)

    for lines in zmat[49:]:

        if "t19" in lines:
            outputf.write("t19= 180.0 S 36 10.0\n")
        elif "t27" in lines:
            outputf.write("t27= %d.0 F\n"%psi)
        else:
            outputf.write(lines)


#############
#MAIN SCRIPT#
#############
inputfiles=glob.glob("gzmat_files/*.gzmat")


if not os.path.exists("submission_files/"):
    os.makedirs("submission_files")


print("Creating  folders and files")
for f in inputfiles:
    reader = open(f,"r").readlines()
    name = f.split("/")[1].split(".gzmat")[0]
    psi_angle = int(name.split("_")[1].split(".")[0]) - 180
    name = "structure_%d" %psi_angle
    directory ="submission_files/%d/" % psi_angle
    os.makedirs(directory)
    ext = ".sh"
    script_name = directory+name+ext
    print("Creating bash script for")
    script_file = open(script_name,"w")
    write_submit(script_file,name)
    ext=".gzmat"
    outputname=directory+name+ext
    outputf=open(outputname,"w")
    print("Creating gzmat for gaussian")
    config_zmat(reader,outputf,psi_angle)
