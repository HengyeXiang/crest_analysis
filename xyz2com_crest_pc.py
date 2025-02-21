import os
import numpy as np
import math
import sys
import matplotlib.pyplot as plt

####################################
# This script is aimed to address the following issue: Assuming we have a .xyz file containing a lot of frames,
# like after conformational sampling using CREST package, we can obtain such a .xyz file containing all conformers 
# within a specified energy threshold. We want to convert each conformer inside that .xyz file into a Gaussian 16 
# input file. Using existing visualization tools like VMD or PyMOL to open that .xyz file, you still have to operate 
# on each frame, in addition to that, those tools are not supported to generate Gaussian 16 input file directly, you
# have to save each frame into a .xyz or .pdb file manually first. This can take significant longer time when the 
# conformation ensenmble become really large. Using general file convertion tools like OpenBabel directly on such .xyz
# file will not generate the gaussian 16 input file for each frame, it will end with a single conformer .xyz file instead.
# Even though you can acheive that with OpenBabel, you still have to manually specify the keywords used in Gaussian 16 calculation
# for each conformer, or write another script to do that. This kind of file generation issue can become quite important to
# optimizing automation workflows.
#
# Herein, I wrote a python script to generate Gaussian 16 input files for each frame (conformer) inside an ensemble .xyz file,
# keywords of Gaussian 16 input files will be provided via a "input.txt" file under the current directory, which I believe is a 
# more flexible way to specify the basic as well as advanced settings we commonly used for Gaussian calculations, compared to 
# specifying a lot of calculation related variables inside this script, since the needs for different people can vary too much,
# if I want to cover them all, it would make the script become too complicated and not easy to use. This approach can be applied 
# to converting each frame (conformer) into any other type file format used in modern Quantum Chemistry calculation packages like 
# Q-Chem or ORCA based on similar idea, as long as we're familiar with what the input file would look like.

# In addition to the file convertion, the script also provided another functionality to calculate the bond distances/bond angles/dihedrals
# for each frame (conformer) based on specified atom indexes (-l for distance, -a for angle and -d for dihedral). It will save these
# geometric related values into a separate .txt file under the generated new folder, as well as showing the scatter plot for them
# in the terminal and saving those plots under current working directory. This functionality is similar to cpptraj module, which is
# commonly used for analyzing MD trajectory results.
#####################################

# function for calculating bond length
def calc_length(array,i,j):
    x = np.array(array[i].split()[1:],dtype=np.float32)
    y = np.array(array[j].split()[1:],dtype=np.float32)
    sum_d = 0
    for v in range(len(x)):
        sum_d += (x[v]-y[v])**2
    distance = math.sqrt(sum_d)
    return round(distance,4)

# function for calculating bond angle
def calc_angle(array,i,j,k):
    atomA = np.array(array[i].split()[1:],dtype=np.float32)
    atomB = np.array(array[j].split()[1:],dtype=np.float32)
    atomC = np.array(array[k].split()[1:],dtype=np.float32)
    AB = atomB - atomA
    BC = atomB - atomC
    length_AB = math.sqrt(sum(v**2 for v in AB))
    length_BC = math.sqrt(sum(v**2 for v in BC))
    cos_angle = np.dot(AB,BC)/(length_AB * length_BC)
    angle = math.acos(cos_angle)
    angle = angle * 180 / math.pi
    return round(angle,4)

# function for calculating dihedral
def calc_dihedral(array,i,j,k,m):
    atomA = np.array(array[i].split()[1:],dtype=np.float32)
    atomB = np.array(array[j].split()[1:],dtype=np.float32)
    atomC = np.array(array[k].split()[1:],dtype=np.float32)
    atomD = np.array(array[m].split()[1:],dtype=np.float32)
    AB = atomB - atomA
    BC = atomC - atomB
    CD = atomD - atomC
    n1 = np.cross(AB,BC)
    n2 = np.cross(BC,CD)
    norm_n1 = math.sqrt(sum(v**2 for v in n1))
    norm_n2 = math.sqrt(sum(v**2 for v in n2))
    cos_dihedral = np.dot(n1,n2)/(norm_n1 * norm_n2)
    dihedral = math.acos(cos_dihedral)
    dihedral = dihedral * 180 / math.pi
    
    n3 = np.cross(n1,n2)
    if np.dot(n3,BC) < 0:
        dihedral = -dihedral
    
    return round(dihedral,4)

def main(path,folder_name):
    with open(path,'r') as f:
        data = f.readlines()
        line_s = str(data[0]).strip() 
        # number of atoms for the molecule
        num_atom = int(line_s.strip('\n')) 
        # number of conformers from CREST
        num_conf = int(len(data)/(num_atom+2)) 
        # threshold of conformers' number we want for file conversion
        trsd = num_conf 

        # check whether user wants to specify g16 input keywords instead of using defaults
        length = len(sys.argv) 
        
        # lists for geometric properties storage
        bd_length_asb = []
        angle_asb = []
        dihedral_asb = []
        
        # indexes used to store bond distances related atoms
        bd_atm1 = 0
        bd_atm2 = 0
        
        # indexes used to store bond angles related atoms
        angle_atm1 = 0
        angle_atm2 = 0
        angle_atm3 = 0
        
        # indexes used to store dihedral related atoms
        dihedral_atm1 = 0
        dihedral_atm2 = 0
        dihedral_atm3 = 0
        dihedral_atm4 = 0
        
        # number of CPUs for g16 calculation
        core = '%nprocshared=24\n' 
        # number of memory for g16 calculation
        mem = '%mem=48GB\n' 
        # whether you need to save checkpoint file
        bool_chk = False 
        
        # list used to store input keywords info
        inpMB = []
        
        # file containing g16 input keyword info
        input_filename = 'input.txt' 
        input_filepath = os.path.join(os.getcwd(),input_filename)
        
        if os.path.isfile(input_filepath):
            # read from existing g16 input keyword file
            with open(input_filepath,'r') as f5:
                info = f5.readlines()
                inpMB = info
        else:
            # default setting of method and basis set without provided keywords
            inpMB.append('# opt freq b3lyp/def2svp em=gd3bj\n')
            inpMB.append('\n')
            inpMB.append('Title\n')
            inpMB.append('\n')
            inpMB.append('0 1\n') 
            
        if (length > 1):
                for k in range(length):
                    # read bond length keyword
                    if(sys.argv[k] == '-l'):
                        bd_atm1 = int(sys.argv[k+1])+1
                        bd_atm2 = int(sys.argv[k+2])+1 
                    # read bond angle keyword
                    elif(sys.argv[k] == '-a'):
                        angle_atm1 = int(sys.argv[k+1])+1
                        angle_atm2 = int(sys.argv[k+2])+1
                        angle_atm3 = int(sys.argv[k+3])+1 
                    # read dihedral keyword
                    elif(sys.argv[k] == '-d'):
                        dihedral_atm1 = int(sys.argv[k+1])+1
                        dihedral_atm2 = int(sys.argv[k+2])+1
                        dihedral_atm3 = int(sys.argv[k+3])+1
                        dihedral_atm4 = int(sys.argv[k+4])+1 
                    # read number of conformers converted keyword
                    elif(sys.argv[k] == '-n'):
                        trsd = int(sys.argv[k+1]) 
                    # read core info keyword
                    elif(sys.argv[k] == '-c'):
                        core = '%nprocshared=' + sys.argv[k+1] + '\n' 
                    # read memory keyword
                    elif(sys.argv[k] == '-m'):
                        mem = '%mem=' + sys.argv[k+1] + 'GB\n' 
                    # read .chk file keyword
                    elif(sys.argv[k] == '-chk'):
                        bool_chk = True 
        
        # check if the saved folder already exists
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        else:
            print(f'Folder {folder_name} already exists.')
            sys.exit()
            
        for i in range(trsd):
            # write bond distances info based on atom pairs if specified
            if(bd_atm1 >=2 and bd_atm2 >=2):
                path_bd_length = './' + folder_name + '/bond_length.txt'
                bd_length_new = calc_length(data,bd_atm1+i*(num_atom+2),bd_atm2+i*(num_atom+2))
                bd_length_asb.append(bd_length_new)
                with open(path_bd_length,'a') as f2:
                    f2.write(str(i+1) + '\t')
                    f2.write(str(bd_length_new))
                    f2.write('\n')
                    
             # write bond angles info based on atom pairs if specified        
            if(angle_atm1 >=2 and angle_atm2 >=2 and angle_atm3 >=2):
                path_angle = './' + folder_name + '/angle.txt'
                angle_new = calc_angle(data,angle_atm1+i*(num_atom+2),angle_atm2+i*(num_atom+2),angle_atm3+i*(num_atom+2))
                angle_asb.append(angle_new)
                with open(path_angle,'a') as f3:
                    f3.write(str(i+1) + '\t')
                    f3.write(str(angle_new))
                    f3.write('\n')
                    
            # write dihedrals info based on atom pairs if specified        
            if(dihedral_atm1 >=2 and dihedral_atm2 >=2 and dihedral_atm3 >=2 and dihedral_atm4 >=2):
                path_dihedral = './' + folder_name + '/dihedral.txt'
                dihedral_new = calc_dihedral(data,dihedral_atm1+i*(num_atom+2),dihedral_atm2+i*(num_atom+2),dihedral_atm3+i*(num_atom+2),dihedral_atm4+i*(num_atom+2))
                dihedral_asb.append(dihedral_new)
                with open(path_dihedral,'a') as f4:
                    f4.write(str(i+1) + '\t')
                    f4.write(str(dihedral_new))
                    f4.write('\n')
            
            # file name generated, assigned index based on original order
            path_save = './' + folder_name + "/crest_conformers_" + str(i+1) + '.com'
            
            # part controlling writing the g16 input files
            with open(path_save,'w+') as f1:
                # write memory keyword
                f1.write(mem) 
                # write CPU keyword
                f1.write(core) 
                
                if(bool_chk == True):
                    # write checkpoint file keyword
                    path_save_chk = os.getcwd() + '/' + folder_name + "/crest_conformers_" + str(i+1) + '.chk'
                    f1.write('%chk=' + path_save_chk + '\n')
                    
                if(len(inpMB) <= 5):
                    # cases only writing basic parts for g16 input files 
                    # usually just first 5 lines 
                    # including calculation keyword, notes, charge and multiplicity settings
                    for k in range(len(inpMB)):
                        f1.write(inpMB[k])
                    for j in range(num_atom):
                        f1.write(data[2+j+i*(num_atom+2)])
                    f1.write('\n')
                    
                elif(len(inpMB) > 5):
                    # cases also writing advanced parts for g16 input files
                    # including settings like constrain, scan or ECP basis sets
                    for k in range(5):
                        f1.write(inpMB[k])
                    for j in range(num_atom):
                        f1.write(data[2+j+i*(num_atom+2)])
                    for m in range(len(inpMB)-5):
                        f1.write(inpMB[m+5])
                    f1.write('\n')

        if(len(bd_length_asb) != 0):
            # draw plot for bond length changes
            x = np.array([i for i in range(len(bd_length_asb))])
            plt.scatter(x,bd_length_asb, color='blue', marker='o')
            plt.title('Bond length changes bewteen ' + str(bd_atm1-1) + ' and ' + str(bd_atm2-1) )
            plt.savefig('Bond_'+str(bd_atm1-1)+'and'+str(bd_atm2-1)+'.png', dpi=1000)
            plt.show()
            
        if(len(angle_asb) != 0):
            # draw plot for bond angle changes
            x = np.array([i for i in range(len(angle_asb))])
            plt.scatter(x,angle_asb, color='blue', marker='o')
            plt.title('Angle changes among ' + str(angle_atm1-1) + ', ' + str(angle_atm2-1) + ' and ' +str(angle_atm3-1))
            plt.savefig('Angle_'+str(angle_atm1-1)+'_'+str(angle_atm2-1)+'and'+str(angle_atm3-1)+'.png', dpi=1000)
            plt.show() 
            
        if(len(dihedral_asb) != 0):
            # draw plot for dihedral changes
            x = np.array([i for i in range(len(dihedral_asb))])
            plt.scatter(x,dihedral_asb, color='blue', marker='o')
            plt.title('Dihedral changes among ' + str(dihedral_atm1-1) + ', ' + str(dihedral_atm2-1) + ', ' + str(dihedral_atm3-1) + ' and ' + str(dihedral_atm4-1))
            plt.savefig('Dihedral_'+str(dihedral_atm1-1)+'_'+str(dihedral_atm2-1)+'_'+str(dihedral_atm3-1)+'and'+str(dihedral_atm4-1)+'.png', dpi=1000)
            plt.show() 
    
    
if __name__ == "__main__":
    # This is usually the default file name generated from CREST run
    path = './crest_conformers.xyz' # Replace with your actual file path for other needs
    folder_name = 'crest_conformers' # folder name containing generated g16 input files
    main(path, folder_name)                                

