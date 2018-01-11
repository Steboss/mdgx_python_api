#Refinement of coordinate.py script for converting gaussian minimized coordinate structures to mdcrd and crd files
# USe:  python gaussian_to_crd.py gaussian_output_file.gout  [future:step of dihedral scan(37)]
# Next development: automatically recognize how manyatoms are present (reading at the beginning of gausisan output) and know then after how many line zmatrix orientation could be
#TODO
#READ n of ATOMS
#READ automatically all the points
# gnerate not only 1d scan but 2d scan
from numpy import *
import math
import re
import os
import sys




def stringsubstitution(coordinates):

    substitute = re.sub("\^","",coordinates)
    substitue2=re.sub('\{', '',substitute)
    substitue3=re.sub('\}','',substitue2)
    substitue4=re.sub('\:','',substitue3)
    substitue5=re.sub('\_','',substitue4)
    substitue6=re.sub(',','',substitue5)
    substitue7=substitue6.replace('\\n','')
    substitue8=substitue7.replace('\'','')


    return substitue8


def mdcrdwriter(listOfLists):

    filename = "ligand_valid_structures.mdcrd"
    coord_file = open(filename,"w")
    coord_file.write("Cpptraj generator")
    counter_coord = 0
    list_index = 0 

    for lista in listOfLists:
        lenght = len(lista)
        last_element = lista[lenght-1]
        print("mdcrd number %d" % list_index)
        list_index+=1
        for num, reading  in enumerate(reader[last_element:last_element+500],1):
            zmatrix = "Z-Matrix orientation:"
            if zmatrix in reading:
                coordinates = str(reader[last_element+num+4 : last_element+num+12])
                string_of_coordinates = stringsubstitution(coordinates)
              
                big_file_index =0

                for coordinate in string_of_coordinates.split()[1:]:

                    if counter_coord <3:
                        counter_coord+=1
                    elif counter_coord == 5:    
                        counter_coord = 0
                        coord = re.sub("\]","", coordinate) 
                        
                        if not big_file_index%10:
                            coord_file.write("\n")
                            big_file_index = 0 
                            if float(coord)<0: 
                                coord_file.write("  %.3f" % float(coord))
                            elif float(coord)>=10:
                                coord_file.write("  %.3f" % float(coord))
                            else:
                                coord_file.write("   %.3f" % float(coord))
                        else:
                            if float(coord)<0: 
                                coord_file.write("  %.3f" % float(coord))
                            elif float(coord)>=10:
                                coord_file.write("  %.3f" % float(coord))
                            else:
                                coord_file.write("   %.3f" % float(coord))
                        big_file_index+=1
                                
                    else:
                        coord = re.sub("\]","", coordinate)

                        counter_coord+=1
                        
                        
                        if not big_file_index%10:
                            coord_file.write("\n")
                            big_file_index = 0 
                            if float(coord)<0: 
                                coord_file.write("  %.3f" % float(coord))
                            elif float(coord)>=10:
                                coord_file.write("  %.3f" % float(coord))
                            else:
                                coord_file.write("   %.3f" % float(coord))
                        else:
                            if float(coord)<0: 
                                coord_file.write("  %.3f" % float(coord))
                            elif float(coord)>=10:
                                coord_file.write("  %.3f" % float(coord))
                            else:
                                coord_file.write("   %.3f" % float(coord))
                        big_file_index+=1
                                                 
        

def crdwriter(listOfLists):


     index_list = 0
     counter_coord = 0
     
     for lista in listOfLists:
         index_list+=1
         #print(index_list)
         lenght = len(lista)
         last_element = lista[lenght-1]

         
         for num, reading  in enumerate(reader[last_element:last_element+500],1):
             zmatrix = "Z-Matrix orientation:"
             if zmatrix in reading:
                 file_name = "amber_inps/ligand" + "_" + str(index_list) +  ".mdcrd"   #for the moment let's call ligand the outputfile
                 
    
                 coord_file = open(file_name, "w")
                 coord_file.write("LIG\n")
                 coord_file.write("     8") #future fund how many atoms in te first line of gaussian
                                  
                 coordinates = str(reader[last_element+num+4 : last_element+num+12])
                 string_of_coordinates = stringsubstitution(coordinates)
              
                 small_file_index =0
                 print("Writing file %s" % file_name)
                 for coordinate in string_of_coordinates.split()[1:]:

                     if counter_coord <3:
                         counter_coord+=1
                     elif counter_coord == 5:    
                         counter_coord = 0
                         coord = re.sub("\]","", coordinate) 
                        
                         if not small_file_index%6:
                             coord_file.write("\n")
                             small_file_index = 0 
                             if float(coord)<0: 
                                 coord_file.write("  %.7f" % float(coord))
                             elif float(coord)>=10:
                                 coord_file.write("  %.7f" % float(coord))
                             else:
                                 coord_file.write("   %.7f" % float(coord))
                         else:
                             if float(coord)<0: 
                                 coord_file.write("  %.7f" % float(coord))
                             elif float(coord)>=10:
                                 coord_file.write("  %.7f" % float(coord))
                             else:
                                 coord_file.write("   %.7f" % float(coord))
                         small_file_index+=1
                                
                     else:
                         coord = re.sub("\]","", coordinate)
 
                         counter_coord+=1
                        
                         
                         if not small_file_index%6:
                             coord_file.write("\n")
                             small_file_index = 0 
                             if float(coord)<0: 
                                 coord_file.write("  %.7f" % float(coord))
                             elif float(coord)>=10:
                                 coord_file.write("  %.7f" % float(coord))
                             else:
                                 coord_file.write("   %.7f" % float(coord))
                         else:
                             if float(coord)<0: 
                                 coord_file.write("  %.7f" % float(coord))
                             elif float(coord)>=10:
                                 coord_file.write("  %.7f" % float(coord))
                             else:
                                 coord_file.write("   %.7f" % float(coord))
                         small_file_index+=1


###########################main###########################

output_crd_folder = "amber_inps"
if not os.path.exists(output_crd_folder):
    os.makedirs(output_crd_folder)


gaussian_file = sys.argv[1]
gaussian = open(gaussian_file, "r")
gaussian.seek(0)
reader = gaussian.readlines()

counter = []
lines_found = [] 

listOfLists = [[] for i in range(37)] #list of lists for every step of the scan




 

print("Finding step of minimization")
for i in range(1,38):
    for num,reading in enumerate(reader,1):
#       
        stringa = " " + str(i) + " out of    37"
        
        if stringa in reading:
            print("Found %s" % stringa)
            
            listOfLists[i-1].append(num)


print("mdcrd file writing")
mdcrdwriter(listOfLists)
print("crd files writing")
crdwriter(listOfLists)            

    

            
