#!/usr/bin/python
# 
# This scripts works only with the prmtop format introduced in Amber7
#


from sys import argv
from sys import exit

# DEFAULT VALUES
get_heavy=False         #get only those torsions involving heavy atoms
get_sidechain=False     #get only those torsions that are not backbone
get_backbone=False      #get only phi,psi and omega angles
select=False		#Select only the angles of a given set of residues

# READING OPTIONS
c=0 #counter
if len(argv[0:])==1:
  print 'Usage: get_tor.py [OPTIONS] amber_topology.top'
  print 'for more info use the option -help'
  exit()
for i in range(1,len(argv[0:])):
  if argv[i]=='-help' :
      print '\nUsage of get_tor.py:'
      print '\nSINOPSIS:'
      print   '           get_tor.py [OPTIONS] amber_topology.top [> input_of_ptraj]'
      print '\nOPTIONS:'
      print   '          -h        Select only torsion angles involving only heavy atoms'
      print   '                    (i.e., not only Hydrogen)'
      print '\n          -b        Select only torsion angles defining the backbone'
      print   '                    conformation of polipeptide molecules'
      print   '                    (Note that -h -b is partially redundant)'
      print '\n          -s        Opposite to -b, select only the torsions that are not involved'
      print   '                    in backbone conformation (i.e., side-chain torsions)'
      print   '                    (Note that -b and -s are incompatible options)'
      print '\n          -sel      SELECTION'
      print   '                    Select only torsions involving a given set of residues, for example:'
      print   '                    -sel 5-20        selects torsions in residues from 5 to 20'
      print   '                    -sel 1,5-20,33   as above, but including two more residues: 1 and 33'
      print '\n          -help     Print this quick help.'
      print '\nEXAMPLES:'
      print   '>> get_tor.py peptide6.top > input_torsions.ptraj'
      print   '>> get_tor.py -h -s -sel 3,4,7-10  peptide4.top > input.ptraj'
      exit()
  if argv[i]=='-sel':
     Selection_arg=argv[i+1]
     select=True
     c=c+2
  if argv[i]=='-b' :
     get_backbone=True
     c=c+1
  if argv[i]=='-s' :
     get_sidechain=True
     c=c+1
  if argv[i]=='-h' :
     get_heavy=True
     c=c+1
  if i==(len(argv[0:])-1):
     topology=open(argv[i],"r")
     c=c+1

# CHECKING OPTIONS
if (len(argv[0:])-1)!=c:
  print '\nERROR1: Check the options'
  exit()
if get_backbone==True and get_sidechain==True:
  print '\nERROR2: Incompatible options (-s) and (-b)'
  print 'If you want to select both sidechain and backbone, do not specify anything'
  exit()

print 'trajin PUT-YOUR-TRAYECTORY-HERE'

#MAKING THE LIST OF SELECTED RESIDUES
temp_list1=[]
temp_list2=[]
SelectedResidues=[]
if (select):
   temp_list1=Selection_arg.split(',')
   for i in range(len(temp_list1)):
     temp_list2=temp_list1[i].split('-')
     if (len(temp_list2)>2):
        print "Bad Argument",temp_list2
        exit()
     elif (len(temp_list2)==2):
        for j in range( int(temp_list2[0]),(int(temp_list2[1])+1) ):
	  SelectedResidues.append(j)
     else:
        SelectedResidues.append(int(temp_list1[i]))
        

# READING THE TOPOLOGY
begin_dih_noH=False
begin_dih_withH=False
begin_atom_name=False
begin_atom_name=False
begin_residue_p=False
atom_name=[]
temp_list=[]
list=[]
res_pointer=[]
    #Reading the atom names
for line in topology:
    temp_list=line.split()
    if line.find('ATOM_NAME') != -1:
        begin_atom_name=True
    if begin_atom_name==True and (line[0]!='%'):
        for i in range(len(temp_list)):
          while len(temp_list[i]) > 4:
            atom_name.append(temp_list[i][0:4])
            temp_list[i]=temp_list[i][4:]
          atom_name.append(temp_list[i])
    if line.find('FLAG CHARGE') !=-1:
       break
    #Reading the residue pointer
for line in topology:
    temp_list=line.split()
    if line.find('RESIDUE_POINTER') != -1:
        begin_residue_p=True
    if begin_residue_p and (line[0]!='%'):
        for i in range(len(temp_list)):
          res_pointer.append(int(temp_list[i]))
    if line.find('BOND_FORCE_CONSTANT') !=-1:
        res_pointer.append(len(atom_name)+1)
        break
    #Reading torsions without hydrogen
for line in topology:
    if line.find('DIHEDRALS_WITHOUT_HYDROGEN') != -1:
         begin_dih_noH=True
    if begin_dih_noH==True and (line[0]!='%'):
       if len(line.split())==10 and (int(line.split()[3])>=0):
           list.append(line.split()[0:4])
       if len(line.split())==10 and (int(line.split()[8])>=0):
           list.append(line.split()[5:9])
       if len(line.split())==5 and (int(line.split()[3])>=0):
           list.append(line.split()[0:4])
    if line.find('EXCLUDED_ATOMS_LIST') !=-1:
       break
    #Reading torsions including hydrogens
topology.seek(0)
if (get_heavy==False):
 for line in topology:
    if line.find('DIHEDRALS_INC_HYDROGEN') != -1:
         begin_dih_withH=True
    if begin_dih_withH==True and (line[0]!='%'):
       if len(line.split())==10 and (int(line.split()[3])>=0):
           list.append(line.split()[0:4])
       if len(line.split())==10 and (int(line.split()[8])>=0):
           list.append(line.split()[5:9])
       if len(line.split())==5 and (int(line.split()[3])>=0):
           list.append(line.split()[0:4])
    if line.find('DIHEDRALS_WITHOUT_HYDROGEN') !=-1:
       break

topology.close()

# GETTING RESNUMER[ATOMNUMBER]
ResNum=[]
for i in range(len(res_pointer)-1):
	for j in range(res_pointer[i],res_pointer[i+1]):
           ResNum.append(i+1)
SelectedAtoms=[]
for i in SelectedResidues :
    for j in range(len(atom_name)):
	if ResNum[j]==i:
          SelectedAtoms.append(j+1)  #Making the list of selected atoms

# MAKING THE FIRST LIST  OF TORSIONS
for i in range(len(list)):           #Making the first list of torsions
	list[i][0]=abs(int(list[i][0]))/3+1
        list[i][1]=abs(int(list[i][1]))/3+1
        list[i][2]=abs(int(list[i][2]))/3+1
        list[i][3]=abs(int(list[i][3]))/3+1

# OBTAINING A NEW LIST WITHOUT REDUNDANCY
new_list=[]
k=0
for i in range(len(list)):
        belong=False
        belongs2back=False
	if (atom_name[list[i][1]-1]=='C' and atom_name[list[i][2]-1]=='N') or \
           (atom_name[list[i][2]-1]=='C' and atom_name[list[i][1]-1]=='N') or \
           (atom_name[list[i][1]-1]=='N' and atom_name[list[i][2]-1]=='CA') or \
           (atom_name[list[i][2]-1]=='N' and atom_name[list[i][1]-1]=='CA') or \
           (atom_name[list[i][1]-1]=='C' and atom_name[list[i][2]-1]=='CA') or \
           (atom_name[list[i][2]-1]=='C' and atom_name[list[i][1]-1]=='CA'):
           belongs2back=True
        if k>=0:
	 for j in range(len(new_list)): #Looking for redundancy in list
	  if (list[i][1]==new_list[j][1] and list[i][2]==new_list[j][2]) or \
             (list[i][1]==new_list[j][2] and list[i][2]==new_list[j][1]):
		belong=True
	if (not belong):                #Obtaining a new list without redundancy
           if (get_backbone) and belongs2back :
	     new_list.append(list[i])
             k+=1
           elif (get_sidechain) and (not belongs2back):
	     new_list.append(list[i])
             k+=1
           elif (not get_backbone) and (not get_sidechain):
             new_list.append(list[i])
             k+=1

# PICKING UP THE SELECTED ATOMS WITH THE OPTIION "-sel"
if select:
   new_list2=[]
   for i in range(len(new_list)):
	belongs1=False
	belongs2=False
	for j in SelectedAtoms:
          if new_list[i][1]==j:
            belongs1=True
          if new_list[i][2]==j:
            belongs2=True
          if belongs1 and belongs2:
            new_list2.append(new_list[i])
            break
   new_list=[]   
   new_list=new_list2   

	        

# PRINTING OUT 
for i in range(len(new_list)):
    	new_list[i][0]='@'+str(new_list[i][0])
    	new_list[i][1]='@'+str(new_list[i][1])
    	new_list[i][2]='@'+str(new_list[i][2])
    	new_list[i][3]='@'+str(new_list[i][3])

for i in range(len(new_list)):
    if (i+1)<10:                    name='d'+'000'+str(i+1)+'.dat'
    if (i+1)>=10 and (i+1)<100:     name='d'+'00'+str(i+1)+'.dat'
    if (i+1)>=100 and (i+1)<1000:   name='d'+'0'+str(i+1)+'.dat'
    if (i+1)>=1000:                 name='d'+str(i+1)+'.dat'
    print 'dihedral ',name[0:5],' ',
    print new_list[i][0],new_list[i][1],new_list[i][2],new_list[i][3], ' out ',name

info=open("atoms_in_tor.info","w")
for i in range(len(new_list)):
    info.write(str(new_list[i][1])[1:]+'  '+str(new_list[i][2])[1:]+'\n')

print 'matrix dist :*@* :*@* out distance_matrix.dat'
