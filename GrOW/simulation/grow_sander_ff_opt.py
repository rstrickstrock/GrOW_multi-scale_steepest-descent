#!/usr/local/python/python-3.6.1/bin/python3
import os
import os.path
import sys
import shutil
import glob
import subprocess
import fileinput

''' REQUIRED ENVIRONMENTAL VARIABLE, DIRECTORIES AND FILES
A. export OUTPATH=      ( e.g. /home/karl/32_CoSMoS/99_Marco/Debug_1/octane_test_1)
B. export working_dir=  (e.g. /home/karl/32_CoSMoS/99_Marco/Debug_1/octane_test_1)
C. export BINDIR=       (e.g. /tmp/test)
D. export TPDUMMY=      (e.g. /home/karl/32_CoSMoS/99_Marco/Debug_1/ExTrM.template.dat)
E. export HELPSCRIPT=   (e.g. /home/karl/32_CoSMoS/99_Marco/Debug_1/replace_placeholders.sh)

1. BINDIR/00_qm_opt/molecule-*-gam.inp.log
2. leaprc.extrm
3. leaprc.extrm.w2p
4. molec.extrm.bcc.mol2
5. Force-field template file (e.g. ExTrM.template.dat)

RUN: grow_sander_ff_opt.py -s1 1.099939 -e1 0.553945 -s2 0.145963 -e2 0.003650
'''

print('\033[94m',"   Executing:", os.path.basename(__file__), '\n', '\x1b[0m')

###########################################################
## Check environment
def OUTPATH_env():
    if 'OUTPATH' in os.environ:
        print('\033[92m','        OUTPATH is set to {}'.format(os.environ['OUTPATH']), '\n', '\x1b[0m')
    else:
        print('\033[92m','        OUTPATH has not been set.', '\x1b[0m', '\n')
        sys.exit()

def BINDIR_env():
    if 'BINDIR' in os.environ:
        print('\033[92m','        BINDIR is set to {}'.format(os.environ['BINDIR']), '\n', '\x1b[0m')
    else:
        print('\033[92m','        BINDIR has not been set.', '\x1b[0m', '\n')
        sys.exit()

def TPDUMMY_env():
    if 'TPDUMMY' in os.environ:
        print('\033[92m','        TPDUMMY is set to {}'.format(os.environ['TPDUMMY']), '\n', '\x1b[0m')
    else:
        print('\033[92m','        TPDUMMY has not been set.', '\x1b[0m', '\n')
        sys.exit()

def HELPSCRIPT_env():
    if 'HELPSCRIPT' in os.environ:
        print('\033[92m','        HELPSCRIPT is set to {}'.format(os.environ['HELPSCRIPT']), '\n', '\x1b[0m')
    else:
        print('\033[92m','        HELPSCRIPT has not been set.', '\x1b[0m', '\n')
        sys.exit()

###########################################################
## Check for other mandatory files
def DOES_FILE_EXIST(MYFILE):
    if not os.path.isfile(MYFILE):
        print('\033[92m','        ' + MYFILE + ' does not exists or is set incorrectly.', '\x1b[0m', '\n')
        sys.exit()

def DOES_DIR_EXIST(DIRECTORY):
    if not os.path.isdir(DIRECTORY):
        print('\033[92m','        The directory ' + DIRECTORY + ' does not exists.', '\x1b[0m', '\n')
        sys.exit()

###########################################################
def relative_energy(completeList = [], *args):
    RELATIVE_E = []
    energylist = []
    for line in completeList:
        #print('line.split(\' \')[1] =',line.split(' ')[1])
        energylist.append(line.split(' ')[1])
    EMIN = energylist[0]
    #EMIN = min(energylist)
    print('First entry in completeList:',completeList[0])
    print('EMIN = ',EMIN)
    for E in range(0,len(energylist)):
        #print('completeList[E].split(\' \')[0] =',completeList[E].split(' ')[0])
        #print('int(energylist[E]) =',float(energylist[E]))
        #print('int(EMIN)) =',float(EMIN))
        RELATIVE_E.append(completeList[E].split(' ')[0] + ' ' + str(float(energylist[E])-float(EMIN)))
    return RELATIVE_E

## Collects final energies of MM minimizations
def GET_AMBER_ENERGY(LOG,moleculeName):
    START = 'FINAL RESULTS'
    END = 'BOND'
    LINES= []
    with open(LOG) as input_data:
        for line in input_data:
            if line.strip() == START:
                break
        # Reads text until the end of the block:
        for line in input_data:
            if line.strip() == END:
                break
            LINES.append(line.split())
    LINES = [x for x in LINES if x]  ## remove empty lists
    #print('LINES: %s' %(LINES))
    AMBER_ENERGIES.append(moleculeName + " " + LINES[1][1])

def create_mm_mol2_coords(INFILE, MOL2):
    START = '@<TRIPOS>ATOM'
    END = '@<TRIPOS>BOND'
    COORD = []
    ATOMLABEL = []
    RESNAME = []
    ATOMTYPES = []
    CHARGES = []
    LINE = []
    ## Extract unique xyz coordinates
    if MOL2.endswith('.mol2'):
        with open(MOL2) as input_data_1:
            for line in input_data_1:
                if line.strip() == START:
                    break
            # Reads text until the end of the block:
            for line in input_data_1:
                if line.strip() == END:
                    break
                COORD.append(line[19:46])

    if MOL2.endswith('.xyz'):
        i=1
        with open(MOL2) as input_data_1:
            for line in input_data_1:
                if i > 2:
                    COORD.append(line[16:68])
                i+=1

    ## Extract proper labels, atom types and charges
    with open(INFILE) as input_data:
        for line in input_data:
            if line.strip() == START:
                break
        for line in input_data:
            if line.strip() == END:
                break
            ATOMLABEL.append(line[0:18])
            RESNAME.append(line[50:65])
            ATOMTYPES.append(line.split()[5])
            CHARGES.append(line.split()[8])

    ## combine them
    for a, b, c, d, e in zip(ATOMLABEL,COORD,ATOMTYPES,RESNAME,CHARGES):
        LINE.append(a + ' ' + b + ' ' + c + ' ' + d + ' ' + e)
    return LINE

def psi2xyz(INFILE,OUTFILE):
    START = 'Final optimized geometry and variables:'
    END = 'Cleaning optimization helper files.'
    GEOM = []
    i = 1
    with open(INFILE) as input_data:
        for line in input_data:
            if line.strip() == START:
                break
        # Reads text until the end of the block:
        for line in input_data:  ## This keeps reading the file
            if (i > 5):          ## Skip 5 lines before recording
                if line.strip() == END:
                    break
                GEOM.append(line)
            #else:
            #    print line
            i += 1

    if GEOM[len(GEOM)-1] == '\n':
        GEOM = GEOM[0:len(GEOM)-1]
    number_of_atoms = len(GEOM)
    #print("number of atoms: %s" %(str(number_of_atoms)))
    filename = os.path.basename(INFILE)
    #print("filename: %s" %(filename))

    f=open(OUTFILE,'w')
    f.write(str(number_of_atoms) + '\n')
    f.write(filename + '\n')
    for coord in GEOM:
        if coord.endswith('\n'):
            f.write(coord)
        else:
            f.write(coord + '\n')
    f.close()


###########################################################
## Check and assign environment
OUTPATH_env()
OUTPATH=os.environ['OUTPATH']

BINDIR_env()
BINDIR=os.environ['BINDIR']

TPDUMMY_env()
TPDUMMY=os.environ['TPDUMMY']

HELPSCRIPT_env()
HELPSCRIPT=os.environ['HELPSCRIPT']

## Check if required directories exist
###### !!! Change this to different directories, because copying, writing, reading 
###### !!! the same files when parallel instances of this script are called, cause trouble 
###### !!! and simulations fail
if not os.path.exists(BINDIR):
    os.makedirs(BINDIR)
if not os.path.exists(BINDIR + '/06_mm_opt/'):
    os.makedirs(BINDIR + '/06_mm_opt/')

BINDIR_QM_MM = os.path.join(os.path.dirname(OUTPATH),'BINDIR')
if not os.path.exists(BINDIR_QM_MM):
    os.makedirs(BINDIR_QM_MM)
if not os.path.exists(os.path.join(BINDIR_QM_MM, '06_mm_opt/')):
    os.makedirs(os.path.join(BINDIR_QM_MM, '06_mm_opt/'))

os.environ['BINDIR_QM_MM'] = BINDIR_QM_MM

DOES_DIR_EXIST(BINDIR + '/00_qm_opt/')
###########################################################
## Variables
DOES_FILE_EXIST(os.environ['MOL2'])
DOES_FILE_EXIST(os.environ['LEAPRC'])
DOES_FILE_EXIST(os.environ['W2P'])

shutil.copy2(os.environ['MOL2'],os.path.join(BINDIR_QM_MM, '06_mm_opt/'))
shutil.copy2(os.environ['LEAPRC'],os.path.join(BINDIR_QM_MM, '06_mm_opt/'))
shutil.copy2(os.environ['W2P'],os.path.join(BINDIR_QM_MM, '06_mm_opt/'))

MOL2_TEMPLATE = os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'),os.path.basename(os.environ['MOL2']))
FF_SOURCE_1 = os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'),os.path.basename(os.environ['LEAPRC']))
FF_SOURCE_2 = os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'),os.path.basename(os.environ['W2P']))
#FF_TEMPLATE_FILE = 'ExTrM.template.dat'

shutil.copy2(os.path.join(os.path.join(BINDIR, '06_mm_opt/'),'ExTrM.Amber.hydrocarbons.dat'),os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'),'ExTrM.Amber.hydrocarbons.dat'))
shutil.copy2(os.path.join(os.path.join(BINDIR, '06_mm_opt/'),'frcmod.extrm.w2p'),os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'),'frcmod.extrm.w2p'))

HEAD = []
TAIL = []
AMBER_ENERGIES = []

###########################################################
## Arguments of system call
args = sys.argv

global SIGMA_1
global EPSILON_1
global SIGMA_2
global EPSILON_2

SIGMA_1 = args[1]
SIGMA_2 = args[2]
EPSILON_1 = args[3]
EPSILON_2 = args[4]

def PARSER_FILE():
    import argparse

    global SIGMA_1
    global EPSILON_1
    global SIGMA_2
    global EPSILON_2

    ## Parser for linecommand
    parser = argparse.ArgumentParser(description="Execute sander using two sets of LJ parameters.")
    parser.add_argument('-s1', dest='SIGMA_1', help="Sigma parameter 1.", required=True, default='1.0', type=float)
    parser.add_argument('-e1', dest='EPSILON_1', help="Epsilon parameter 1.", required=True, default='0.01', type=float)
    parser.add_argument('-s2', dest='SIGMA_2', help="Sigma parameter 2.", required=True, default='1.1', type=float)
    parser.add_argument('-e2', dest='EPSILON_2', help="Epsilon parameter 2. (defaults to Gamess", required=True, default='0.05', type=float)
    args = parser.parse_args()

    SIGMA_1 = args.SIGMA_1
    EPSILON_1 = args.EPSILON_1
    SIGMA_2 = args.SIGMA_2
    EPSILON_2 = args.EPSILON_2

###########################################################
## Initialize variables and parse command line

#PARSER_FILE()

## Call helper script
SUBDIR_PATH=OUTPATH + "/subdir"

os.environ["working_dir"] = OUTPATH

if not os.path.exists(SUBDIR_PATH):
    os.makedirs(SUBDIR_PATH)

if not os.path.exists(os.path.join(BINDIR_QM_MM, '06_mm_opt/')):
    os.makedirs(os.path.join(BINDIR_QM_MM, '06_mm_opt/'))

subprocess.call(HELPSCRIPT + ' %s %s %s %s' % (SIGMA_1,SIGMA_2,EPSILON_1,EPSILON_2), shell=True)

## Copy over the template force-field file
shutil.copy2(TPDUMMY,BINDIR_QM_MM + "/06_mm_opt/")

#QMLOG = glob.glob(BINDIR + '/00_qm_opt/molecule-*-gam.inp.log')    ## original was looking for -psi.inp.log
QMLOG = glob.glob(BINDIR + '/00_qm_opt/molecule-*-psi.inp.log')    ## original was looking for -psi.inp.log
#QMLOG = glob.glob(BINDIR + '/00_qm_opt/molecule-00[1-5]-gam.inp.log')    ## For testing
QMLOG=sorted(QMLOG)

###########################################################
## Gather the AMBER proper head and end of mol2 file
with open(BINDIR_QM_MM + '/06_mm_opt/molec.extrm.bcc.mol2') as file:
    HEAD = [next(file) for x in range(6)]
    for line in file:
        if '@<TRIPOS>BOND' in line:
            for line in file:
                TAIL.append(line)

HEAD = [s.rstrip() for s in HEAD]
TAIL = [s.rstrip() for s in TAIL]

## Amber files needed for referencing the force force.
## Replace location of w2p force-field supplementary parameters.
shutil.copy2(os.environ['W2P'],os.path.join(BINDIR_QM_MM, '06_mm_opt/'))
shutil.copy2(os.environ['LEAPRC'],os.path.join(BINDIR_QM_MM, '06_mm_opt/'))
for line in fileinput.input([os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'), os.path.basename(os.environ['W2P']))], inplace=True):
    print(line.replace('SOURCEDIR', BINDIR_QM_MM + '/06_mm_opt'), end='')
for line in fileinput.input([os.path.join(os.path.join(BINDIR_QM_MM, '06_mm_opt/'), os.path.basename(os.environ['LEAPRC']))], inplace=True):
    print(line.replace('SOURCEDIR', BINDIR_QM_MM + '/06_mm_opt'), end='')

###########################################################
## 1. Create a mol2 file from the QM log files such that AMBER's tleap can read it.
## 2. Create leap input and run tleap
## 3. Execute Sander minimization
## 4. Convert optimized geometry rst file to pdb for easy viewing
## 5. Grab raw energy out of the output file.

with open(BINDIR_QM_MM + '/06_mm_opt/min.in', 'w') as f:
    f.write('Constraint Minimization\n')
    f.write('&cntrl\n')
    f.write('imin=1, dielc=1,ntb=0,\n')
    f.write('maxcyc=20000, cut=40.0,\n')
    f.write('ntc=1, ntf=1,\n')
    f.write('drms=0.01,nmropt=0\n')
    f.write('&end\n')


for LOG in QMLOG:
    LOGNAME = LOG.split('/')
    BASENAME = (LOGNAME[-1]).split('.')
    ## For use with GAMESS QM log files
#    BABELGAMMES = ['babel', ' -gamout', LOG, BASENAME[0] + '.xyz']
#    BABELGAMMES = ['babel', ' -gamout', LOG, BASENAME[0] + '.mol2']
    psi2xyz(LOG, os.path.join(OUTPATH, BASENAME[0] + '.xyz'))
#    BABELXYZ = ['babel', XYZ, BASENAME[0] + '.mol2']

    ## Create a mol2 file from a GAMESS log file
    FNULL = open(os.devnull, 'w')
#    subprocess.call(BABELGAMMES, stdout=FNULL, stderr=subprocess.STDOUT)

    ## Collect unique coordinates for each structure
    MOL2_COORDS = []
#    MOL2_COORDS = create_mm_mol2_coords(MOL2_TEMPLATE, BASENAME[0] + '.mol2')
    MOL2_COORDS = create_mm_mol2_coords(MOL2_TEMPLATE, os.path.join(OUTPATH, BASENAME[0] + '.xyz'))

    ## Put it all together into a single file
#    with open(BASENAME[0] + '.mol2', 'w') as f:
    with open(os.path.join(OUTPATH,BASENAME[0] + '.mol2'),'w') as f:
        for item in HEAD:
            f.write("%s\n" % item)
        for item in MOL2_COORDS:
            f.write("%s\n" % item)
        f.write("@<TRIPOS>BOND\n")
        for item in TAIL:
            f.write("%s\n" % item)
#    shutil.copy2(BASENAME[0] + '.mol2',BINDIR + '/06_mm_opt')
    shutil.copy2(os.path.join(OUTPATH,BASENAME[0] + '.mol2'),BINDIR_QM_MM + '/06_mm_opt')

    ## Create each leap.in
    with open(OUTPATH + BASENAME[0] + '.leap.in', 'w') as f:
        f.write('logfile ' + OUTPATH + BASENAME[0] + '.leap.log\n')
        f.write('source ' + FF_SOURCE_1 +'\n')
        f.write('source ' + FF_SOURCE_2 +'\n')
        f.write('verbosity 2\n')
        f.write('a = loadmol2 ' + os.path.join(BINDIR_QM_MM, '06_mm_opt/') + BASENAME[0] + '.mol2\n')
        #f.write('saveamberparm a ' + BINDIR + '/06_mm_opt/' + BASENAME[0] + '.leap.top ' + BINDIR + '/06_mm_opt/' + BASENAME[0] + '.leap.crd\n')
        f.write('saveamberparm a ' + OUTPATH + BASENAME[0] + '.leap.top ' + OUTPATH + BASENAME[0] + '.leap.crd\n')
        f.write('savepdb a ' + os.path.join(BINDIR_QM_MM, '06_mm_opt/') + BASENAME[0] + '.leap.pdb\n')
        f.write('quit\n')

    ## Run tleap to get filename.top and filename.crd --> used as input for MM minimization or MD
    TLEAP = ['tleap', '-s -f', OUTPATH + BASENAME[0] + '.leap.in']
    subprocess.call(TLEAP, stdout=FNULL, stderr=subprocess.STDOUT)
    DOES_FILE_EXIST(OUTPATH + BASENAME[0] + '.leap.top')

    ## Run minimizations
    SANDER = ['sander', '-O', '-i', os.path.join(BINDIR_QM_MM, '06_mm_opt/') + 'min.in', '-o', OUTPATH + BASENAME[0] + '.min.out', '-p', OUTPATH + BASENAME[0] + '.leap.top', '-c', OUTPATH + BASENAME[0] + '.leap.crd', '-r', OUTPATH + BASENAME[0] + '.min.rst', '-ref', OUTPATH + BASENAME[0] + '.leap.crd']
    #SANDER = ['sander', '-O', '-i', BINDIR + '/06_mm_opt/' + 'min.in', '-o', BINDIR + '/06_mm_opt/' + BASENAME[0] + '.min.out', '-p', BINDIR + '/06_mm_opt/' + BASENAME[0] + '.leap.top', '-c', BINDIR + '/06_mm_opt/' + BASENAME[0] + '.leap.crd', '-r', BINDIR + '/06_mm_opt/' + BASENAME[0] + '.min.rst', '-ref', BINDIR + '/06_mm_opt/' + BASENAME[0] + '.leap.crd']
    subprocess.call(SANDER, stdout=FNULL, stderr=subprocess.STDOUT)
    #DOES_FILE_EXIST(BINDIR + '/06_mm_opt/' + BASENAME[0] + '.min.out')
    DOES_FILE_EXIST(OUTPATH + BASENAME[0] + '.min.out')
    
    ## create pdb from MM minimization restart file (i.e. final optimized structure).
    AMBPDB = ['ambpdb', '-p', OUTPATH + BASENAME[0] + '.leap.top', '-c', OUTPATH + BASENAME[0] + '.min.rst']
    #process = subprocess.Popen(AMBPDB, stdin=subprocess.PIPE, stdout=open(BINDIR + '/06_mm_opt/' + BASENAME[0] + '.min.rst.pdb', 'w'), stderr=subprocess.PIPE)
    subprocess.call(AMBPDB, stdout=open(OUTPATH + BASENAME[0] + '.min.rst.pdb', 'w'), stderr=FNULL)

    ## Obtain raw energy
    #GET_AMBER_ENERGY(BINDIR + '/06_mm_opt/' + BASENAME[0] + '.min.out')
    GET_AMBER_ENERGY(OUTPATH + BASENAME[0] + '.min.out',BASENAME[0])
###########################################################
## write out raw energy data


#AMBER_ENERGIES = [float(i) for i in AMBER_ENERGIES]
#AMBER_ENERGIES.sort()

#with open(BINDIR + '/06_mm_opt/Energy.raw.extrm.txt', 'w') as f:
with open(OUTPATH + '/Energy.raw.extrm.txt', 'w') as f:
    for item in AMBER_ENERGIES:
        f.write('%s\n' % item)


RELATIVE_E = relative_energy(AMBER_ENERGIES)

## write out relative energy data
#with open(BINDIR + '/06_mm_opt/Energy.rel.extrm.txt', 'w') as f:
with open(OUTPATH + '/Energy.rel.extrm.txt', 'w') as f:
    for item in RELATIVE_E:
        f.write('%s\n' % item)

# HARDCODED!!!
f = open('/scratch/rstric2s/current_sim/octane_molecule_target_names.txt', 'r')
test_names = f.readlines()
f.close()

f = open(OUTPATH + '/Energy.rel.extrm.txt', 'r')
rel_e = f.readlines()
f.close()

if not len(test_names) == len(rel_e):
    print('len(test_names) = %s and len(rel_e) = %s are not equal!' %(len(test_names), len(rel_e)))
    sys.exit()

properties=[]
for line in range(0,len(test_names)):
    if not test_names[line].split('.')[0] == rel_e[line].split(' ')[0].split('.')[0]:
        print('molecule names %s and %s are not equal, might be a wrong ordering.' %(test_names[line].split('.')[0], rel_e[line].split(' ')[0].split('.')[0]))
        sys.exit()
    properties.append(str(rel_e[line].split(' ')[1]))
        
with open(OUTPATH + '/properties.txt','w') as f:
    for item in properties:
        if not item.endswith('\n'):
            f.write('%s\n' % item)
        else:
            f.write('%s' % item)
#os.system("cp %s %s" %(BINDIR + '/06_mm_opt/Energy.rel.extrm.txt', OUTPATH + '/properties.txt'))
#os.system("cp %s %s" %(OUTPATH + '/Energy.rel.extrm.txt', OUTPATH + '/properties.txt'))
#print(OUTPATH + '/properties.txt')
os.chdir(OUTPATH)
os.system("cp properties.txt ..")

if os.path.isfile(OUTPATH + '/properties.txt'):
    os.system("echo %s > %s" %("Simulation successfully terminated.",OUTPATH + '/terminated.txt'))
else:
    print("grow_sander_ff_opt.py: %s not found. Application exits." %(OUTPATH + '/properties.txt'))
    sys.exit()
    
## Exit
print('\033[94m', "   Exiting:", os.path.basename(__file__), '\n','\x1b[0m')
sys.exit()
