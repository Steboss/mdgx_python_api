#June 2016 Stefano Bosisio
#Script to extract z-matrix or crd files from gaussian output
#gout comes from scan
#Desirable output@ z-matrix with optimization keywords for gaussian jobs
#and coordinate files for Sire

import os,sys,re

def merge_gout(first,second):
    first = open(sys.argv[1]).readlines()
    second = open(sys.argv[2]).readlines()
    outputfile = open("output.gout","w")
    print("Reading first output...")
    for f in first:
        outputfile.write(f)
    print("Reading second output...")
    for f in second:
        outputfile.write(f)
    print("create file output.gout")
    outputfile.close()

def create_template(reader,start,end):
    #Creation of a template for zmatrix recovery

    zmatrix_template = open("z_mat_template.dat","w")
    for zmat in reader[start:end+1]: #this is from teh gaussian output, so re-formatting is needed
        splitter=zmat.split()
        #print(splitter)
        #print(len(splitter))
        if len(splitter)<3:  #should be only 1 atom
            zmatrix_template.write("    %s\n" % splitter[0])
        elif len(splitter) <4: #should containt one bond
            zmatrix_template.write("    %s    %s      %s\n" %(splitter[0],splitter[1],splitter[2]))
        elif len(splitter)< 6: #should have bond + angle
            zmatrix_template.write("    %s    %s      %s    %s      %s\n" %(splitter[0],splitter[1],splitter[2],splitter[3],splitter[4]))

        else : #should have torsional too
            atom=splitter[0]
            number =splitter[1]
            bond = splitter[2]
            bondnumb = splitter[3]
            angle = splitter[4]
            anglenumb = splitter[5]
            torsional = splitter[6]

            if int(number)>9:
                numb="   %s" %number
            else:
                numb="    %s" %number


            if (int((re.findall(r'(\d+)',bond)[0]))) > 9:
                #print("Number is %d and is 9" %(int((re.findall(r'(\d+)',bond)[0]))))
                bondvar="     %s" %bond
            else:
                bondvar="      %s" %bond

            if int(bondnumb)>9:
                bondvarnumb="   %s" %bondnumb
            else:
                bondvarnumb="    %s" %bondnumb

            if (int((re.findall(r'(\d+)',angle)[0]))) > 9:
                anglevar="     %s" %angle
            else:
                anglevar="      %s" %angle

            if int(anglenumb)>9:
                anglevarnumb="   %s" %anglenumb
            else:
                anglevarnumb="    %s" %anglenumb

            if (int((re.findall(r'(\d+)',torsional)[0]))) > 9:
                torsvar="     %s" %torsional
            else:
                torsvar="      %s"%torsional

            content="    %s%s%s%s%s%s%s\n" %(atom,numb,bondvar,bondvarnumb,anglevar,anglevarnumb,torsvar)
            zmatrix_template.write(content)

    zmatrix_template.close()



def find_zmatrix(reader,stationary_index,natoms,template,variable):

    last_idx = len(stationary_index)
    end_reader=len(reader)

    for i in range(0,last_idx):
        #From here we can do a function and we pass the stationary index
        #Here we want the first two Gradgrad indexes, where the structure is containted
        counter=0
        #print("GradGradGrad lines search")

        starting_point = stationary_index[i]
        if i==last_idx-1:
            ending_point= end_reader
        else:
            ending_point=stationary_index[i+1]

        extract=reader[starting_point:ending_point]
        print("Finding Z-matrix")
        #the first occurrence is wat we need
        zmat_idx = extract.index('                         Z-Matrix orientation:                         \n')

        crd_start = zmat_idx + 5
        crd_end = zmat_idx + 5 + natoms

        zmatrix = extract[crd_start:crd_end]
        #first write the zmatix, so we can know the values of variable/degree angle
        phi =  gzmatwriter(extract,natoms,variable,template)

        crdwriter(zmatrix,natoms,phi)



def crdwriter(zmatrix,natoms,phi):

    #Crd file

    outcrdname = "crd/structure_%s.crd"  %phi
    outputcrd = open(outcrdname,"w")
    outputcrd.write("LIG\n")
    outputcrd.write("    %d\n" %natoms)

    #Now write the coordinates rst7 file

    position = 0
    counter = 0
    for f in zmatrix:
        coordx = float(f.split()[3])
        coordy = float(f.split()[4])
        coordz = float(f.split()[5])

        if coordx<0 or coordx>10:
            space="  "
            crdX = "%.7f" % coordx
            cX = space+crdX
        else:
            space="   "
            crdX = "%.7f" % coordx
            cX = space+crdX

        if coordy<0 or coordy>10:
            space="  "
            crdY = "%.7f" % coordy
            cY = space+crdY
        else:
            space="   "
            crdY = "%.7f" % coordy
            cY = space+crdY

        if coordz<0 or coordz>10:
            space="  "
            crdZ = "%.7f" % coordz
            cZ = space+crdZ
        else:
            space="   "
            crdZ = "%.7f" % coordz
            cZ = space+crdZ

        if counter ==1 :
            outputcrd.write("%s%s%s\n" %(cX,cY,cZ))
            counter=0
        else:
            outputcrd.write("%s%s%s" %(cX,cY,cZ))
            counter+=1

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

def gzmatwriter(reader,natoms,variable,template):

    #Necessary step to understand at what anle we are
    #Sometimes Gaussian do not peform the minimization and so structures are skipped

    opt_start = reader.index('                       !   Optimized Parameters   !\n') +5
    opt_end = reader.index(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n')-1
    dict_variable = {}
    for param in reader[opt_start:opt_end]:
        variable = param.split()[1]
        value = param.split()[2]
        dict_variable[variable]=(value)

    for key in list(dict_variable.keys()):
        if key==variable:
            phi = (dict_variable[key])



    gzmat_deg_folder = "gzmat/%s" %phi
    if not os.path.exists(gzmat_deg_folder):#
        os.makedirs(gzmat_deg_folder)
    #here we have   extract_gzmat/ANGLE_0/structure_10.gzmat ...
    outgzmatname ="%s/structure_%s.gzmat" %(gzmat_deg_folder,phi)
    #bash script to be save inthe same folde of gzmat
    bashname = "%s/structure_%s.sh" %(gzmat_deg_folder,phi)
    outputbash = open(bashname,"w")
    structure_name = "structure_%s" %phi
    write_submit(outputbash,structure_name)
    ######################################
    outputgzmat = open(outgzmatname,"w")
    outputgzmat.write("%chk=structure.chk\n")
    outputgzmat.write("NProcShared=8\n")
    outputgzmat.write("#HF/6-31G* SCF=tight opt=Z-Matrix Freq\n")
    outputgzmat.write("\n")
    outputgzmat.write("Optimization of phi %s\n" %phi )
    outputgzmat.write("\n")
    outputgzmat.write("0   1\n")

    #Zmatrix optimized values are immediately after the stataionary point
    opt_start = reader.index('                       !   Optimized Parameters   !\n') +5
    opt_end = reader.index(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n')-1
    #print(reader[opt_start])
    for temp in template:
        outputgzmat.write(temp)
    outputgzmat.write("Variables\n")

    for param in reader[opt_start:opt_end]:
        variable = param.split()[1]
        value = param.split()[2]
        if float(value)<0 and float(value)>-10:
            outputgzmat.write("%s=  %s\n"%(variable,value))
        elif float(value)<-10:
            outputgzmat.write("%s= %s\n"%(variable,value))
        elif float(value)>9 and float(value)<99:
            outputgzmat.write("%s=  %s\n"%(variable,value))
        elif float(value)>99:
            outputgzmat.write("%s= %s\n"%(variable,value))
        else:
            outputgzmat.write("%s=   %s\n" %(variable,value))

    return phi

    #parmasplit
#        outputgzmat.write(param)


#############
#MAIN SCRIPT#
#############


if len(sys.argv)==3:
    print("Merging the two output")
    merge_gout(sys.argv[1],sys.argv[2])
    inputfile = open("output.gout","r")
else:
    print("One output given")
    inputfile = open(sys.argv[1],"r")


reader = inputfile.readlines()
counter = 0
#Let's find out how many ATOMS we have
for i,f in enumerate(reader,0):
    if "Charge" in f:
        if counter==0: #first charge should be found
            start = i+1
            counter+=1
        else:
            continue
    elif "Variables" in f: #then variables
        if counter==1:
            end = i-1
            break  #thus break
print(start,end)
number_atoms = end-start+1
print("The number of atoms is: %d" % number_atoms)
######################################################
#Now create a zmatrix TEMPLATE
#print(reader[start:end])
create_template(reader,start,end)

template=open("z_mat_template.dat","r").readlines()
###################################################
#Now find the variables that is scanned
var_start = reader.index( '       Variables:\n')
var_end = reader.index(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n')
print("let's find out the variable to be optimized")
for f in reader[var_start:var_end]:
    if "Scan" in f:
        variable = f.split()[0]  #variable which is scanned
    else:
        continue
print("The variable is %s is it correct?" % variable)
####################################################
#Now find out Stationary POINT positions
stationary_index = []
print("Stationary point:")
for i,f in enumerate(reader,0):
    if "Stationary point" in f:
        stationary_index.append(i)
    else:
        continue
print(stationary_index)
#########################################################

#Create a saving directroy
outputdir = "crd"
if not os.path.exists(outputdir):#
    os.makedirs(outputdir)
outputdir = "gzmat"
if not os.path.exists(outputdir):#
    os.makedirs(outputdir)

find_zmatrix(reader,stationary_index,number_atoms,template,variable)
