These scripts can become quite useful after you finished CREST runs and wants to do some analysis based on the generated conformer ensemble.

Main functionalities: 
1) file conversion of each conformer into gaussian 16 input file;
2) calculate bond distances/bond angles/dihedrals for each conformer, save values into .txt file, generate and save corresponding scatter plots;

   For using the script, try the following commands:
   python3 xyz2com_crest_pc.py -l 1 2 (calculate distance between atom 1 and 2)
   python3 xyz2com_crest_pc.py -a 1 2 3 (calculate angle between atom 1,2 and 3)
   python3 xyz2com_crest_pc.py -d 1 2 3 4 (calculate dihedral among atom 1,2,3 and 4)
   The keywords can be combined together (orders are not important), like
   python3 xyz2com_crest_pc.py -d 1 2 3 4 -l 1 2 -a 1 2 3
   
4) calculate RMSD values for each conformer related to the first one (lowest-energy one from CREST), save those data into .txt file.

Example input files and output files are provided as well.
