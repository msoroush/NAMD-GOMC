import os
import datetime
import subprocess
import pandas as pd
import scipy as sp
import numpy as np

from warnings import warn

# MD/MC simulations (NAMD/GOMC) for the NPT, NVT, GEMC, and GCMC ensembles
#*************************************************
# Changable variables (start)
#*************************************************
# simulations to performs.
#Note: at this time the energy is not checked between the ending simulation and the restarted simulation
Total_cycles_NAMD_GOMC_sims = 1000 # 400000 # enter integer: # simulations NAMD + GOMC sims = 2 * Total_cycles_NAMD_GOMC_sims
Starting_at_cycle_NAMD_GOMC_sims = 1    # enter integer: enter 0 for intial simualtion start, other for restart

#  GPU/CPU selection, Simulation type (ensemble), cores in each box,
Use_CPU_or_GPU =  'GPU'   # enter string : 'CPU' or 'GPU', currently only CPU available,  you need to provide the NAMD exectutable GPU code as well
simulation_type = 'NVT'  # enter string : 'GEMC', 'GCMC', 'NPT', 'NVT' (only GEMC-NVT available currently: 'GEMC' = GEMC-NVT)
only_use_box_0_for_NAMD_for_GEMC = True    # bool, True or False.  True = NAMD runs only box 0 for GEMC, False = NAMD runs box 0 and 1 for GEMC
No_core_box_0 =  6  # CPU cores: this is the ONLY place to enter CPU cores for 'GCMC', 'NPT', 'NVT', and  'GEMC' and only_use_box_0_for_NAMD_for_GEMC = True
No_core_box_1 =  0  # CPU cores:  is always ZERO for 'GCMC', 'NPT', 'NVT'.  Only use when 'GEMC' and only_use_box_0_for_NAMD_for_GEMC = True

#Physical/Chemical properties Properties
simulation_Temp_K = 500 # float: GOMC and NAMD units of temperature are in K
simulation_Pressure_bar =  None # 1.01325 bar = 1 atm # float: GOMC and NAMD units of pressure are in bar (must have number even if unused)
# Currently, a total of only 4 chemicals/residues can be used for GCMC simulations
GCMC_only_ChemPot_Type_1_Residue = None  # string: GCMC ONLY!!! -- GOMC residue name of the molecule, must have Type_1 for GCMC
GCMC_only_ChemPot_Type_1_mu_K = None  # float: GCMC ONLY!!! -- GOMC units of mu (chemical potential ) are in units of K, must have Type_1 for GCMC
GCMC_only_ChemPot_Type_2_Residue =  None  # string: GCMC ONLY!!! -- GOMC residue name of the molecule or None if not used
GCMC_only_ChemPot_Type_2_mu_K = None  # float: GCMC ONLY!!! -- GOMC units of mu (chemical potential ) are in units of K or None if not used
GCMC_only_ChemPot_Type_3_Residue = None  # string: GCMC ONLY!!! -- GOMC residue name of the molecule or None if not used
GCMC_only_ChemPot_Type_3_mu_K = None  # float: GCMC ONLY!!! -- GOMC units of mu (chemical potential ) are in units of K or None if not used
GCMC_only_ChemPot_Type_4_Residue = None  # string: GCMC ONLY!!! -- GOMC residue name of the molecule or None if not used
GCMC_only_ChemPot_Type_4_mu_K = None  # float: GCMC ONLY!!! -- GOMC units of mu (chemical potential ) are in units of K or None if not used

GOMC_Run_Steps = 5000 #4000  # needs to be 10 minimum for now, NEEDS TO BE THE SAME AS THE PREVIOUS SIMULATION, IF RESTARTED!
NAMD_Run_Steps = 5000 #800   # needs to be 10 minimum for now, NEEDS TO BE THE SAME AS THE PREVIOUS SIMULATION, IF RESTARTED!
NAMD_Minimize_mult_scalar = 0   # interger, NAMD_minimize steps = NAMD_Run_Steps * NAMD_Minimize_mult_scalar

# BOX 0 set_x_dim, set_y_dim, set_z_dim  if set to None this reads the data from the PDB files CRYST1 line
set_x_dim_box_0 =  None # float or interger, or None if provided in PDB, units in Angstrom.
set_y_dim_box_0  = None  # float or interger, or None if provided in PDB units in Angstrom.
set_z_dim_box_0  = None  # float or interger, or None if provided in PDB units in Angstrom.

# BOX 1 set_x_dim, set_y_dim, set_z_dim  if set to None this reads the data from the PDB files CRYST1 line
# for GEMC and GCMC ONLY !!!!!!!!
set_x_dim_box_1 =  None # float or interger, or None if provided in PDB units in Angstrom.
set_y_dim_box_1  = None # float or interger, or None if provided in PDB units in Angstrom.
set_z_dim_box_1  = None # float or interger, or None if provided in PDB units in Angstrom.

# read_angle_alpha, read_angle_beta, read_angle_gamma  if set to None this reads the data from the PDB files CRYST1 line
# ENTER ONLY  None or 90 , as only othoganal box is currently accepted!!!!!!!!!!!
set_angle_alpha = None # float or interger, angle in Degrees
set_angle_beta = None # float or interger, angle in Degrees.
set_angle_gamma = None # float or interger, angle in Degrees.

Starting_FF_file_GOMC = 'required_data/equilb_box_500K/GOMC_pore_water_rigid_FF.inp'  # from main folder
Starting_FF_file_NAMD = 'required_data/equilb_box_500K/NAMD_pore_water_FF.inp'  # from main folder
# box 0 data
Starting_PDB_box_0_file = '/required_data/equilb_box_500K/Equilb_500K_20ns_pore_3x3x2.8nm_1-layer.pdb' # from main folder
Starting_PSF_box_0_file = '/required_data/equilb_box_500K/Equilb_500K_20ns_pore_3x3x2.8nm_1-layer.psf' # from main folder

# box 1 data
Starting_PDB_box_1_file = 'required_data/equilb_box_500K/reservior_box.pdb' # from main folder
Starting_PSF_box_1_file = 'required_data/equilb_box_500K/reservior_box.psf' # from main folder

# MD/MC engine file names and paths, if the paths are required
NAMD_executable_file = str(os.getcwd()) + '/required_data/bin/NAMD212/namd2' + '_' + str(Use_CPU_or_GPU)
GOMC_executable_file = str(os.getcwd()) + '/required_data/bin/GOMC_'+ str(Use_CPU_or_GPU) + '_'  + str(simulation_type)

allowable_error_fraction_vdw = 1 * 10**(-3)
allowable_error__fraction_electro = 1 * 10**(-3)
max_absolute_allowable_kcal_fraction_electro = 0.5
#*************************************************
# Changable variables (end)
#*************************************************




#*************************************************
# Potential changable variables in the future (start)
#*************************************************

GOMC_console_BLKavg_Hist_Steps = int( GOMC_Run_Steps  )
GOMC_RST_Coor_CKpoint_Steps = int( GOMC_Run_Steps  )

if GOMC_Run_Steps/10 > 500:
    GOMC_Hist_sample_Steps  = int( 500 )
elif GOMC_Run_Steps/10 <= 500:
    GOMC_Hist_sample_Steps = int(GOMC_Run_Steps/10)

NAMD_RST_DCD_XST_Steps = int( NAMD_Run_Steps   )
NAMD_console_BLKavg_E_and_P_Steps = int( NAMD_Run_Steps   )

use_GOMC_checkpoint = 'true' # select 'true' or 'false', needs to be false for now as does not print seperate PSF's


#*************************************************
# Potential changable variables in the future (end)
#*************************************************

#*************************************************
# check for existing and create NAMD and GOMC folders (start)
#*************************************************
if os.path.isdir('NAMD') or os.path.isdir('GOMC'):
    warn('INFORMATION: if the system fails to start (with errors) from the beginning of a simulation, '
         'you may need to delete the main GOMC and NAMD folders.  '
         'The failure to start/restart may be caused by the last simulation not finishing correctly.'
         )
    warn('INFORMATION: If the system fails to restart a previous run (with errors), '
         'you may need to delete the last subfolders under the main '
         'NAMD and GOMC (i.e., folders NAMD = 00000000_a or GOMC = 00000001). '
         'The failure to start/restart may be caused by the last simulation not finishing properly.'
         )

path_NAMD_runs = "NAMD"
os.makedirs(path_NAMD_runs, exist_ok=True)

path_GOMC_runs = "GOMC"
os.makedirs(path_GOMC_runs, exist_ok=True)
#*************************************************
# create NAMD and GOMC folders (end)
#*************************************************

#*************************************************
# create NAMD and GOMC config file templates locations (start)
#*************************************************
path_NAMD_template = 'required_data/config_files/NAMD.conf'

path_GOMC_template = 'required_data/config_files/GOMC_'+str(simulation_type) +'.conf'
#*************************************************
# create NAMD and GOMC config file templates locations (end)
#*************************************************

#check the pressure values and set to atomospheric if not used
if simulation_type in ['NPT'] and (isinstance(simulation_Pressure_bar, float) == False \
                                   and isinstance(simulation_Pressure_bar, int) == False):
    raise Exception("The simulation pressure needs to be set for the NPT simulation type (int or float). \n")
if simulation_type in ['NPT'] and (
        isinstance(simulation_Pressure_bar, float) == True or isinstance(simulation_Pressure_bar, int) == True):
    if simulation_Pressure_bar < 0:
        raise Exception("The simulation pressure needs to be set to a positive or zero [(int of float) >=0 bar] "
                        "value for the NPT simulation type. \n")
if simulation_type in ['NVT', 'GEMC', 'GCMC'] :
    #set pressure to atomospheric as it needs a number and it is not used
    simulation_Pressure_bar = 1.01325



K_to_kcal_mol = 1.98720425864083 * 10**(-3)

Python_file_directory = os.getcwd()

Log_Template_file = open(str(Python_file_directory) + "/NAMD_GOMC_started_at_cycle_No_" \
                         + str(Starting_at_cycle_NAMD_GOMC_sims)+ ".log", 'w')

start_time = datetime.datetime.today()
write_log_data = "*************************************************\n" + \
                 'date and time (start) = ' + str(start_time) + " \n" \
                 "*************************************************" + ' \n'
Log_Template_file.write(str(write_log_data))
print(str(write_log_data))

write_log_data = "*************************************************\n" + \
                 'NAMD_executable_file = ' +str(NAMD_executable_file) + " \n" \
                 "*************************************************" + ' \n'
Log_Template_file.write(str(write_log_data))
print(str(write_log_data))

write_log_data = "*************************************************\n" + \
                 'GOMC_executable_file = ' +str(GOMC_executable_file) + " \n" \
                 "*************************************************" + ' \n'
Log_Template_file.write(str(write_log_data))
print(str(write_log_data))


NAMD_Minimize_Steps = int(int(NAMD_Run_Steps) * int(NAMD_Minimize_mult_scalar))
if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
    Total_No_cores = No_core_box_0 + No_core_box_1

    if No_core_box_1 == 0:
        write_log_data = "*************************************************\n" \
                         + "No_core_box_0 was set to " + str(No_core_box_0) + ",\n" \
                         + "No_core_box_0 = " + str(No_core_box_0) + " \n" \
                         + 'WARNING: the number of CPU cores listed for box 1 is zero (0), and should be an ' \
                         + 'interger > 0, or the NAMD simulation for box 1 will not run and the python script ' \
                         + 'will crash \n' \
                         + "No_core_box_1 was set to " + str(No_core_box_1) + ",\n" \
                         + "No_core_box_1 = " + str(No_core_box_1) + " \n" \
                         + "*************************************************" + " \n"
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

    else:
        write_log_data = "*************************************************\n" \
                         + "No_core_box_0 was set to " + str(No_core_box_0) + ",\n" \
                         + "No_core_box_0 = " + str(No_core_box_0) + " \n" \
                         + "No_core_box_1 was set to " + str(No_core_box_1) + ",\n" \
                         + "No_core_box_1 = " + str(No_core_box_1) + " \n" \
                         + "*************************************************" + " \n"
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

else:
    if No_core_box_1 != 0:
        write_log_data = "*************************************************\n"  \
                         + "No_core_box_0 was set to " + str(No_core_box_0) + ",\n" \
                         + "No_core_box_0 = " + str(No_core_box_0) + " .\n" \
                         + 'WARNING: the number of CPU cores listed for box 1 are not being used, \n' \
                         + "No_core_box_1 was set to " + str(No_core_box_1) + ",\n" \
                         + "No_core_box_1 = " + str(No_core_box_1) + "\n" \
                         + "*************************************************" + " \n"
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))
        No_core_box_1 = 0

    Total_No_cores = No_core_box_0

#Check total cores vs per box cores
if isinstance(No_core_box_0, int) == False:
    write_log_data = "*************************************************\n" + \
                     'Enter No_core_box_0 as an interger, No_core_box_0 = '+ str(No_core_box_0) + " \n"\
                     "*************************************************" + ' \n'
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))
    if No_core_box_0 ==0 :
        write_log_data = "*************************************************\n" + \
                         'Enter No_core_box_0 as a non-zero number, No_core_box_0 = ' + str(No_core_box_0) + " \n" \
                         "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

if isinstance(No_core_box_1, int) == False:
    write_log_data = "*************************************************\n" + \
                     'Enter No_core_box_1 as an interger, No_core_box_1 = '+ str(No_core_box_1) + " \n"\
                     "*************************************************" + ' \n'
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))
    if No_core_box_1 ==0 :
        write_log_data = "*************************************************\n" + \
                         'Enter No_core_box_1 as a non-zero number, No_core_box_1 = ' + str(No_core_box_1) + " \n" \
                         "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

def Calc_folder_zeros(Run_Number):
    total_run_No_digits = 10
    add_No_zeros_at_start_Run_No = int(total_run_No_digits - len(str(Run_Number)))
    add_Number_of_zeros_at_start_Run_No_str = ""
    for z in range(0, add_No_zeros_at_start_Run_No):
        add_Number_of_zeros_at_start_Run_No_str += '0'
    return add_Number_of_zeros_at_start_Run_No_str



# run the GOMC and NAMD simulations

Total_sims_NAMD_GOMC = int(2 * Total_cycles_NAMD_GOMC_sims)
Starting_sims_NAMD_GOMC = int(2 * Starting_at_cycle_NAMD_GOMC_sims)


default_NAMD_E_titles = ['ETITLE:', 'TS', 'BOND', 'ANGLE', 'DIHED', 'IMPRP', 'ELECT', 'VDW', 'BOUNDARY',
                        'MISC', 'KINETIC', 'TOTAL', 'TEMP', 'POTENTIAL', 'TOTAL3', 'TEMPAVG', 'PRESSURE',
                        'GPRESSURE', 'VOLUME', 'PRESSAVG', 'GPRESSAVG']


def get_NAMD_run_0_PME_DIM(box_number):
    # box number can be box 0 or 1
    if  isinstance(box_number,int) == False and  (box_number != 0 or box_number != 1):
        warn( "Enter an interger of 0 or 1  for box_number in the get_NAMD_run_0_PME_DIM function \n" )

    NAMD_first_Run_No = int(0)
    add_zeros_for_box_X_Run_0 = Calc_folder_zeros(NAMD_first_Run_No)

    if box_number == 0:
        NAMD_box_X_Run_0_dir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                               str(add_zeros_for_box_X_Run_0) + str(NAMD_first_Run_No) + "_a"

    elif box_number == 1:
        NAMD_box_X_Run_0_dir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                               str(add_zeros_for_box_X_Run_0) + str(NAMD_first_Run_No) + "_b"

    try:
        read_NAMD_box_X_Run_0_log_file = open(NAMD_box_X_Run_0_dir + "/" \
                                              + "out.dat", 'r').readlines()

        for i, line in enumerate(read_NAMD_box_X_Run_0_log_file):
            split_line = line.split()
            if len(split_line) >= 7:
                if split_line[0] == 'Info:' and split_line[1] == 'PME' and \
                        split_line[2] == 'GRID' and split_line[3] == 'DIMENSIONS':
                    NAMD_X_PME_GRID_DIM = int(split_line[4])
                    NAMD_Y_PME_GRID_DIM = int(split_line[5])
                    NAMD_Z_PME_GRID_DIM = int(split_line[6])

    except:
        NAMD_X_PME_GRID_DIM = None
        NAMD_Y_PME_GRID_DIM = None
        NAMD_Z_PME_GRID_DIM = None

    return NAMD_X_PME_GRID_DIM,  NAMD_Y_PME_GRID_DIM, NAMD_Z_PME_GRID_DIM, NAMD_box_X_Run_0_dir


def get_NAMD_run_0_FFT_filename(box_number):
    # box number can be box 0 or 1
    if isinstance(box_number, int) == False and (box_number != 0 or box_number != 1):
        warn("Enter an interger of 0 or 1  for box_number in the get_NAMD_run_0_PME_DIM function \n")

    NAMD_first_Run_No = int(0)
    add_zeros_for_box_X_Run_0 = Calc_folder_zeros(NAMD_first_Run_No)

    if box_number == 0:
        NAMD_box_X_Run_0_dir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                               str(add_zeros_for_box_X_Run_0) + str(NAMD_first_Run_No) + "_a"

    elif box_number == 1:
        NAMD_box_X_Run_0_dir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                               str(add_zeros_for_box_X_Run_0) + str(NAMD_first_Run_No) + "_b"
    try:
        NAMD_box_X_Run_0_file_list = sorted(os.listdir(NAMD_box_X_Run_0_dir))
    except:
        NAMD_box_X_Run_0_file_list = None

    try:
        for q in range(0, len(NAMD_box_X_Run_0_file_list)):
            NAMD_box_X_Run_0_FFT_1st_9_char_filename = str(NAMD_box_X_Run_0_file_list[q])[0:9]
            if NAMD_box_X_Run_0_FFT_1st_9_char_filename == 'FFTW_NAMD':
                NAMD_box_X_Run_0_FFT_NAMD_filename = NAMD_box_X_Run_0_file_list[q]

        return NAMD_box_X_Run_0_FFT_NAMD_filename

    except:
        write_log_data = "*************************************************\n" + \
                         'The NAMD FFT file not deteted from Run 0 in Box ' + str(box_number) + " \n" \
                         "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))






def delete_NAMD_run_0_FFT_file(box_number):
    # box number can be box 0 or 1
    if  isinstance(box_number,int) == False and  (box_number != 0 or box_number != 1):
        warn( "Enter an interger of 0 or 1  for box_number in the get_NAMD_run_0_PME_DIM function \n" )

    try:
        NAMD_box_X_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number)

        # delete the FFT file in the run 0 box
        rm_NAMD_FFT_box_0_command = "rm " + str(NAMD_box_X_Run_0_dir) + "/" \
                            + str(NAMD_box_X_Run_0_FFT_NAMD_filename)
        exec_rm_NAMD_FFT_box_0_command = subprocess.Popen(rm_NAMD_FFT_box_0_command, shell=True, stderr=subprocess.STDOUT)

        write_log_data = "*************************************************\n" + \
                         'The NAMD FFT file was deleted from Run 0 in Box ' + str(box_number) + " \n" \
                         "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

    except:
        write_log_data = "*************************************************\n" + \
                         'The NAMD FFT file not deteted from Run 0 in Box ' + str(box_number) + " \n" \
                         "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

    return


def check_for_PDB_dims_and_override(DIM_axis, Run_No, read_dim, set_dim = None, only_on_Run_No = 0):
    #read_dim = int or float
    #set_dim  =int or float set by user
    # DIM_axis =  string x y or z
    # only_on_Run_No = int; only overiden on the selected run number

    # check user  override y-dimensions
    if Run_No == only_on_Run_No:
        if read_dim == None:
            if set_dim != None and (isinstance(set_dim, float) or isinstance(set_dim, int)):
                used_dim = set_dim
            else:
                write_log_data = "WARNING: The user defined " + str(DIM_axis) + "-dimension is None " \
                                 + "or not an integer or a float, and the PDB file has no dimension information " \
                                 + ' \n'
                Log_Template_file.write(str(write_log_data))
                warn(str(write_log_data))

        elif read_dim != None and set_dim != read_dim and (isinstance(set_dim, float) or isinstance(set_dim, int)) :
            write_log_data = "WARNING: The user defined " + str(DIM_axis) + "-dimension is different " \
                             + "than the one read from the starting PDB file " \
                             + str(DIM_axis) + "-dim_PDB = " + str(read_dim) + " " + str(DIM_axis) \
                             + "-dim_user_set = " + str(set_dim) \
                             + ". The code is setting the user defined " + str(DIM_axis) + "-dimension.  " \
                             + ' \n'
            Log_Template_file.write(str(write_log_data))
            warn(str(write_log_data))
            used_dim = set_dim

        else:
            used_dim = read_dim

    else :
        used_dim = read_dim

    return used_dim

def write_NAMD_conf_file(Python_file_directory, path_NAMD_template, path_NAMD_runs, GOMC_newdir,
                         Run_No, box_number, NAMD_Run_Steps, NAMD_Minimize_Steps, NAMD_RST_DCD_XST_Steps,
                         NAMD_console_BLKavg_E_and_P_Steps, simulation_Temp_K, simulation_Pressure_bar,
                         Starting_PDB_box_X_file, Starting_PSF_box_X_file,
                         NAMD_X_PME_GRID_DIM, NAMD_Y_PME_GRID_DIM, NAMD_Z_PME_GRID_DIM,
                         set_x_dim = None, set_y_dim= None, set_z_dim= None,
                         set_angle_alpha = 90, set_angle_beta = 90, set_angle_gamma = 90,
                         FFT_add_NAMD_Ang_to_box_dim = 0):

    # The box angles are not currenatly available and are default set to 90 degress.
    # they can be added later, but for now all simulation but be orthoganol boxes

    # box number can be box 0 or 1
    if  isinstance(box_number,int) == False and  (box_number != 0 or box_number != 1):
        warn( "Enter an interger of 0 or 1  for box_number in the get_NAMD_run_0_PME_DIM function \n" )

    # get NAMD box_x directory
    add_zeros_at_start_Run_No_str = Calc_folder_zeros(Run_No)
    if box_number == 0:
        NAMD_box_X_newdir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                            str(add_zeros_at_start_Run_No_str) + str(Run_No) + "_a"
    elif box_number == 1:
        NAMD_box_X_newdir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                            str(add_zeros_at_start_Run_No_str) + str(Run_No) + "_b"
    os.makedirs(NAMD_box_X_newdir, exist_ok=True)
    generate_NAMD_file = open(NAMD_box_X_newdir + "/" + 'in.conf', 'w')

    NAMD_Template_file = open(str(Python_file_directory) + "/" + path_NAMD_template, 'r')
    NAMD_Template_data = NAMD_Template_file.read()
    NAMD_Template_file.close()

    if Run_No != 0:

        write_log_data = "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

        new_NAMD_data = NAMD_Template_data.replace("parameter_file",
                                                   str(Python_file_directory) + "/" + Starting_FF_file_NAMD)
        new_NAMD_data = new_NAMD_data.replace("pdb_box_file",
                                              str(GOMC_newdir) + "/" + "Output_data_BOX_" + str(box_number) \
                                              + "_restart.pdb")
        new_NAMD_data = new_NAMD_data.replace("psf_box_file",
                                              str(GOMC_newdir) + "/" + "Output_data_BOX_" + str(box_number) \
                                              + "_restart.psf")

        new_NAMD_data = new_NAMD_data.replace("coor_file",
                                              str(GOMC_newdir) + "/" + "Output_data_BOX_" \
                                              + str(box_number)+ "_restart.coor")
        new_NAMD_data = new_NAMD_data.replace("xsc_file",
                                              str(GOMC_newdir) + "/" + "Output_data_BOX_"\
                                              + str(box_number)+ "_restart.xsc")
        new_NAMD_data = new_NAMD_data.replace("Bool_restart", str("true"))

        # Read the angles from the intial/starting PDB file
        read_PDB_file = open(str(GOMC_newdir) + "/" + "Output_data_BOX_"
                             + str(box_number)+ "_restart.pdb", 'r').readlines()

    else:

        new_NAMD_data = NAMD_Template_data.replace("parameter_file",
                                                   str(Python_file_directory) + "/" + Starting_FF_file_NAMD)
        new_NAMD_data = new_NAMD_data.replace("pdb_box_file",
                                              str(Python_file_directory) + "/" + Starting_PDB_box_X_file)
        new_NAMD_data = new_NAMD_data.replace("psf_box_file",
                                              str(Python_file_directory) + "/" + Starting_PSF_box_X_file)
        new_NAMD_data = new_NAMD_data.replace("coor_file", str("NA"))
        new_NAMD_data = new_NAMD_data.replace("xsc_file", str("NA"))
        new_NAMD_data = new_NAMD_data.replace("Bool_restart", str("false"))

        read_PDB_file = open(str(Python_file_directory) + "/" + Starting_PDB_box_X_file, 'r').readlines()

    read_x_dim = None
    read_y_dim = None
    read_z_dim = None
    read_angle_alpha = None
    read_angle_beta = None
    read_angle_gamma = None
    # Read the angles from the intial/starting PDB file
    for i, line in enumerate(read_PDB_file):
        if 'CRYST1' in line:
            read_x_dim = read_PDB_file[i].split()[1:2]
            read_x_dim = float(read_x_dim[0])
            read_y_dim = read_PDB_file[i].split()[2:3]
            read_y_dim = float(read_y_dim[0])
            read_z_dim = read_PDB_file[i].split()[3:4]
            read_z_dim = float(read_z_dim[0])

            read_angle_alpha = read_PDB_file[i].split()[4:5]
            read_angle_alpha = float(read_angle_alpha[0])
            read_angle_beta = read_PDB_file[i].split()[5:6]
            read_angle_beta = float(read_angle_beta[0])
            read_angle_gamma = read_PDB_file[i].split()[6:7]
            read_angle_gamma = float(read_angle_gamma[0])

    #check user  override x-dimensions
    used_x_dim = check_for_PDB_dims_and_override("x", Run_No, read_x_dim, set_dim =set_x_dim, only_on_Run_No = 0)

    # check user  override y-dimensions
    used_y_dim = check_for_PDB_dims_and_override("y", Run_No, read_y_dim, set_dim =set_y_dim, only_on_Run_No = 0)

    # check user  override z-dimensions
    used_z_dim = check_for_PDB_dims_and_override( "z", Run_No, read_z_dim, set_dim =set_z_dim, only_on_Run_No = 0)

    # check box for othoganality
    # check the read_angle_alpha
    if read_angle_alpha != None and read_angle_alpha != 90 :
        if Run_No == 0:
            write_log_data = " WARNING: The alpha angle is not 90 degress as read from the " \
                             + "starting PDB file, " \
                             + "read_angle_alpha_PDB = " + str(read_angle_alpha) \
                             + ". Only alpha angle of 90 degress or None is currently allowed.  "  + ' \n'
            Log_Template_file.write(str(write_log_data))
            warn(str(write_log_data))

    if set_angle_alpha != None and set_angle_alpha != 90 :
        if Run_No == 0:
            write_log_data = " WARNING: The alpha angle is not 90 degress as set by the user," \
                             + "set_angle_alpha = " + str(set_angle_alpha) \
                             + ". Only alpha angle of 90 degress or None is currently allowed.  " + ' \n'
        Log_Template_file.write(str(write_log_data))
        warn(str(write_log_data))

    # check the read_angle_beta
    if read_angle_beta != None and read_angle_beta != 90:
        if Run_No == 0:
            write_log_data = " WARNING: The beta angle is not 90 degress as read from the " \
                             + "starting PDB file, " \
                             + "read_angle_beta_PDB = " + str(read_angle_beta) \
                             + ". Only beta angle of 90 degress or None is currently allowed.  " + ' \n'
            Log_Template_file.write(str(write_log_data))
            warn(str(write_log_data))

    if set_angle_beta != None and set_angle_beta != 90:
        if Run_No == 0:
            write_log_data = " WARNING: The beta angle is not 90 degress as set by the user " \
                             + "set_angle_beta = " + str(set_angle_beta) \
                             + ". Only beta angle of 90 degress or None is currently allowed.  " + ' \n'
        Log_Template_file.write(str(write_log_data))
        warn(str(write_log_data))

    # check the read_angle_gamma
    if read_angle_gamma != None and read_angle_gamma != 90:
        if Run_No == 0:
            write_log_data = " WARNING: The gamma angle is not 90 degress as read from the " \
                             + "starting PDB file " \
                             + "read_angle_gamma_PDB = " + str(read_angle_gamma) \
                             + ". Only gamma angle of 90 degress or None is currently allowed.  " + ' \n'
            Log_Template_file.write(str(write_log_data))
            warn(str(write_log_data))

    if set_angle_gamma != None and set_angle_gamma != 90:
        if Run_No == 0:
            write_log_data = " WARNING: The gamma angle is not 90 degress as set by the user " \
                             + "set_angle_gamma = " + str(set_angle_gamma) \
                             + ". Only gamma angle of 90 degress or None is currently allowed.  " + ' \n'
        Log_Template_file.write(str(write_log_data))
        warn(str(write_log_data))

    new_NAMD_data = new_NAMD_data.replace("x_dim_box", str(used_x_dim))
    new_NAMD_data = new_NAMD_data.replace("y_dim_box", str(used_y_dim))
    new_NAMD_data = new_NAMD_data.replace("z_dim_box", str(used_z_dim))
    new_NAMD_data = new_NAMD_data.replace("x_origin_box", str(0))
    new_NAMD_data = new_NAMD_data.replace("y_origin_box", str(0))
    new_NAMD_data = new_NAMD_data.replace("z_origin_box", str(0))

    new_NAMD_data = new_NAMD_data.replace("NAMD_Run_Steps", str((int(NAMD_Run_Steps))))
    new_NAMD_data = new_NAMD_data.replace("NAMD_Minimize", str(int(NAMD_Minimize_Steps)))
    new_NAMD_data = new_NAMD_data.replace("NAMD_RST_DCD_XST_Steps", str((int(NAMD_RST_DCD_XST_Steps))))
    new_NAMD_data = new_NAMD_data.replace("NAMD_console_BLKavg_E_and_P_Steps",
                                          str(int(NAMD_console_BLKavg_E_and_P_Steps)))

    new_NAMD_data = new_NAMD_data.replace("current_step", str((int(0))))
    new_NAMD_data = new_NAMD_data.replace("System_temp_set", str(simulation_Temp_K))
    new_NAMD_data = new_NAMD_data.replace("System_press_set", str(simulation_Pressure_bar))

    if Run_No != 0:
        new_NAMD_data = new_NAMD_data.replace("X_PME_GRID_DIM", str(NAMD_X_PME_GRID_DIM))
        new_NAMD_data = new_NAMD_data.replace("Y_PME_GRID_DIM", str(NAMD_Y_PME_GRID_DIM))
        new_NAMD_data = new_NAMD_data.replace("Z_PME_GRID_DIM", str(NAMD_Z_PME_GRID_DIM))

    else:
        # add x number times more point to the PME grid for "GEMC", "NPT".
        # scalar_dim_mult = 1.3 --> allows for 2x change in volume,
        # allowable volume change = scalar_dim_mult^3 = 1.3^3
        if simulation_type in ["GEMC", "NPT"]:
            scalar_dim_mult = 1.3
            used_and_scaled_NAMD_PME_x_dim = (used_x_dim + FFT_add_NAMD_Ang_to_box_dim)  * scalar_dim_mult
            used_and_scaled_NAMD_PME_y_dim = (used_y_dim + FFT_add_NAMD_Ang_to_box_dim)  * scalar_dim_mult
            used_and_scaled_NAMD_PME_z_dim = (used_z_dim + FFT_add_NAMD_Ang_to_box_dim)  * scalar_dim_mult
        else:
            scalar_dim_mult = 1
            used_and_scaled_NAMD_PME_x_dim = (used_x_dim + FFT_add_NAMD_Ang_to_box_dim) * scalar_dim_mult
            used_and_scaled_NAMD_PME_y_dim = (used_y_dim + FFT_add_NAMD_Ang_to_box_dim) * scalar_dim_mult
            used_and_scaled_NAMD_PME_z_dim = (used_z_dim + FFT_add_NAMD_Ang_to_box_dim) * scalar_dim_mult

        new_NAMD_data = new_NAMD_data.replace("X_PME_GRID_DIM", str(int(used_and_scaled_NAMD_PME_x_dim)))
        new_NAMD_data = new_NAMD_data.replace("Y_PME_GRID_DIM", str(int(used_and_scaled_NAMD_PME_y_dim)))
        new_NAMD_data = new_NAMD_data.replace("Z_PME_GRID_DIM", str(int(used_and_scaled_NAMD_PME_z_dim)))

    generate_NAMD_file.write(new_NAMD_data)
    generate_NAMD_file.close()

    write_log_data = "NAMD simulation data for simulation number " + str(Run_No) + " in box " + str(box_number) \
                                              + " is completed" + ' \n'
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))

    return NAMD_box_X_newdir




def get_NAMD_energy_data(read_NAMD_box_X_Energy_file, E_default_NAMD_titles):
    Get_E_titles = True
    E_values_NAMD_List = []
    for i, line in enumerate(read_NAMD_box_X_Energy_file):
        if line.startswith('ETITLE:') == True and Get_E_titles == True:
            E_titles_NAMD_iteration = read_NAMD_box_X_Energy_file[i]
            E_titles_NAMD_iteration = E_titles_NAMD_iteration.split()
            E_titles_NAMD_iteration = E_titles_NAMD_iteration[:]
            if Get_E_titles == True:
                Get_E_titles = False

        if line.startswith('ENERGY:') == True:
            E_values_NAMD_iteration = read_NAMD_box_X_Energy_file[i]
            E_values_NAMD_iteration = E_values_NAMD_iteration.split()
            E_values_NAMD_List.append(E_values_NAMD_iteration)

    try:
        E_titles_NAMD_iteration
    except:
        E_titles_NAMD_iteration = E_default_NAMD_titles

    NAMD_energy_data_box_X_df = pd.DataFrame(data=E_values_NAMD_List, columns=E_titles_NAMD_iteration)

    # extract energy data from box 0
    NAMD_E_electro_box_X = NAMD_energy_data_box_X_df.loc[:, 'ELECT']
    NAMD_E_electro_box_X_initial_value = float(NAMD_E_electro_box_X.values.tolist()[0])
    NAMD_E_electro_box_X_final_value = float(NAMD_E_electro_box_X.values.tolist()[-1])
    NAMD_E_potential_box_X = NAMD_energy_data_box_X_df.loc[:, 'POTENTIAL']
    NAMD_E_potential_box_X_initial_value = float(NAMD_E_potential_box_X.values.tolist()[0])
    NAMD_E_potential_box_X_final_value = float(NAMD_E_potential_box_X.values.tolist()[-1])
    NAMD_E_vdw_box_X = NAMD_energy_data_box_X_df.loc[:, 'VDW']
    NAMD_E_vdw_box_X_initial_value = float(NAMD_E_vdw_box_X.values.tolist()[0])
    NAMD_E_vdw_box_X_final_value = float(NAMD_E_vdw_box_X.values.tolist()[-1])

    return NAMD_E_electro_box_X, NAMD_E_electro_box_X_initial_value, NAMD_E_electro_box_X_final_value,\
           NAMD_E_potential_box_X, NAMD_E_potential_box_X_initial_value, NAMD_E_potential_box_X_final_value, \
           NAMD_E_vdw_box_X, NAMD_E_vdw_box_X_initial_value, NAMD_E_vdw_box_X_final_value



def compare_NAMD_GOMC_energies(E_vdw_box_0_final_value, E_vdw_box_0_initial_value,
                               E_electro_box_0_final_value, E_electro_box_0_initial_value,
                               Run_No, box_number):
    try:
        error_fract_in_vdw_box_0 = np.abs((E_vdw_box_0_final_value - E_vdw_box_0_initial_value) / \
                                          E_vdw_box_0_final_value)
    except:
        if E_vdw_box_0_final_value == 0 and E_vdw_box_0_initial_value == 0:
            error_fract_in_vdw_box_0 = 0
        else:
            error_fract_in_vdw_box_0 = 'NA'

    if error_fract_in_vdw_box_0 is not 'NA' and error_fract_in_vdw_box_0 <= allowable_error_fraction_vdw:
        write_log_data = "PASSED:  Box " + str(box_number) \
                         + " VDW error fraction between the check between the last point in run " \
                         + str(int(Run_No - 1)) + " and the first point in run " + str(int(Run_No)) \
                         + ", error fraction =  " + str(error_fract_in_vdw_box_0) + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))
    else:
        write_log_data = "FAILED:  Box " + str(box_number) + " VDW error fraction between the last point in run " + \
                         str(int(Run_No - 1)) + " and the first point in run " + str(
            int(Run_No)) + ", error fraction =  " \
                         + str(error_fract_in_vdw_box_0) + ' \n'
        Log_Template_file.write('WARNING: ' + str(write_log_data))
        warn(write_log_data)

    # calc error in Electrostatics box 0
    try:
        error_fract_in_electro_box_0 = np.abs((E_electro_box_0_final_value - E_electro_box_0_initial_value) / \
                                              E_electro_box_0_final_value)
    except:
        if E_electro_box_0_final_value == 0 and E_electro_box_0_initial_value == 0:
            error_fract_in_electro_box_0 = 0
        else:
            error_fract_in_vdw_box_1 = 'NA'

    abs_diff__in_electro_box_0 = np.abs(E_electro_box_0_final_value - E_electro_box_0_initial_value)

    if error_fract_in_electro_box_0 is not 'NA' and error_fract_in_electro_box_0 <= allowable_error__fraction_electro:
        write_log_data = "PASSED: Box " + str(box_number) \
                         + " The Electrostatic error fraction between  the last point in run " + \
                         str(int(Run_No - 1)) + " and the first point in run " + str(int(Run_No)) \
                         + ", error fraction =  " + str(error_fract_in_electro_box_0) + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

    elif abs_diff__in_electro_box_0 <= max_absolute_allowable_kcal_fraction_electro:
        write_log_data = "PASSED: Box " + str(box_number) \
                         + " The Electrostatic error fraction between  the last point in run " + \
                         str(int(Run_No - 1)) + " and the first point in run " + str(int(Run_No)) + \
                         ", absolute difference is =  " + str(abs_diff__in_electro_box_0) + " kcal/mol.  " + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))
    else:
        write_log_data = "FAILED: Box " + str(box_number) \
                         + " Electrostatic error fraction between the last point in run " + \
                         str(int(Run_No - 1)) + " and the first point in run " + str(int(Run_No)) \
                         + ", error fraction =  " + str(error_fract_in_electro_box_0) \
                         + " or the absolute difference is =  " \
                         + str(abs_diff__in_electro_box_0) + " kcal/mol.  " + ' \n'
        Log_Template_file.write('WARNING: ' + str(write_log_data))
        warn(str(write_log_data))


def write_GOMC_conf_file(Python_file_directory, path_GOMC_runs, Run_No, GOMC_Run_Steps,
                         GOMC_RST_Coor_CKpoint_Steps, GOMC_console_BLKavg_Hist_Steps, GOMC_Hist_sample_Steps,
                         simulation_Temp_K, simulation_Pressure_bar,
                         Starting_PDB_box_0_file, Starting_PDB_box_1_file):

    # Create the GOMC configuration file
    GOMC_Template_file = open(str(Python_file_directory) + "/" + path_GOMC_template, 'r')
    GOMC_Template_data = GOMC_Template_file.read()
    GOMC_Template_file.close()
    GOMC_newdir = str(Python_file_directory) + "/" + path_GOMC_runs + "/" \
                  + str(add_zeros_at_start_Run_No_str) + str(Run_No)
    os.makedirs(GOMC_newdir, exist_ok=True)

    write_log_data = "*************************************************" + ' \n'
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))

    generate_GOMC_file = open(GOMC_newdir + "/" + 'in.conf', 'w')

    new_GOMC_data = GOMC_Template_data.replace("parameter_file",
                                               str(Python_file_directory) + "/" + Starting_FF_file_GOMC)

    new_GOMC_data = new_GOMC_data.replace("coor_box_0_file",
                                          NAMD_box_0_newdir + '/namdOut.coor')
    new_GOMC_data = new_GOMC_data.replace("xsc_box_0_file",
                                          NAMD_box_0_newdir + '/namdOut.xsc')




    if previous_GOMC_dir =='NA':
        new_GOMC_data = new_GOMC_data.replace("Restart_Checkpoint_file", 'false')
        new_GOMC_data = new_GOMC_data.replace("pdb_file_box_0_file",
                                              str(Python_file_directory) + "/" + Starting_PDB_box_0_file)
        new_GOMC_data = new_GOMC_data.replace("psf_file_box_0_file",
                                              str(Python_file_directory) + "/" + Starting_PSF_box_0_file)
    else:
        new_GOMC_data = new_GOMC_data.replace("pdb_file_box_0_file",
                                              str(previous_GOMC_dir) + "/" + "Output_data_BOX_0_restart.pdb")
        new_GOMC_data = new_GOMC_data.replace("psf_file_box_0_file",
                                              str(previous_GOMC_dir) + "/" + "Output_data_BOX_0_restart.psf")


    # Read the angles from the intial/starting PDB file
    NAMD_xsc_box_0_file = NAMD_box_0_newdir + "/namdOut.xsc"
    read_NAMD_xsc_box_0_file = open(NAMD_xsc_box_0_file, 'r').readlines()

    read_x_dim_box_0 = read_NAMD_xsc_box_0_file[-1].split()[1:2]
    read_x_dim_box_0 = float(read_x_dim_box_0[0])
    read_y_dim_box_0 = read_NAMD_xsc_box_0_file[-1].split()[5:6]
    read_y_dim_box_0 = float(read_y_dim_box_0[0])
    read_z_dim_box_0 = read_NAMD_xsc_box_0_file[-1].split()[9:10]
    read_z_dim_box_0 = float(read_z_dim_box_0[0])

    read_x_dim_origin_box_0 = read_NAMD_xsc_box_0_file[-1].split()[10:11]
    read_x_dim_origin_box_0 = float(read_x_dim_origin_box_0[0])
    read_y_dim_origin_box_0 = read_NAMD_xsc_box_0_file[-1].split()[11:12]
    read_y_dim_origin_box_0 = float(read_y_dim_origin_box_0[0])
    read_z_dim_origin_box_0 = read_NAMD_xsc_box_0_file[-1].split()[12:13]
    read_z_dim_origin_box_0 = float(read_z_dim_origin_box_0[0])

    new_GOMC_data = new_GOMC_data.replace("x_dim_box_0", str(read_x_dim_box_0))
    new_GOMC_data = new_GOMC_data.replace("y_dim_box_0", str(read_y_dim_box_0))
    new_GOMC_data = new_GOMC_data.replace("z_dim_box_0", str(read_z_dim_box_0))


    if simulation_type in ["GEMC", "GCMC"]:
        readlines_GOMC_Template_file = open(str(Python_file_directory) + "/" + path_GOMC_template, 'r').readlines()
        if simulation_type in ["GCMC"] or (simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == True):
            if previous_GOMC_dir == "NA":
                for i, line in enumerate(readlines_GOMC_Template_file):
                    split_line = line.split()
                    if line.startswith('binCoordinates') == True or line.startswith('extendedSystem') == True:
                        if split_line[0] == 'binCoordinates' and split_line[1] == '1':
                            new_GOMC_data = new_GOMC_data.replace(line, '')
                        elif split_line[0] == 'extendedSystem' and split_line[1] == '1':
                            new_GOMC_data = new_GOMC_data.replace(line, '')

                read_PDB_file = open(str(Python_file_directory) + "/" + Starting_PDB_box_1_file, 'r').readlines()

                read_x_dim_box_1 = None
                read_y_dim_box_1  = None
                read_z_dim_box_1= None
                read_angle_alpha = None
                read_angle_beta = None
                read_angle_gamma = None
                # Read the angles from the intial/starting PDB file
                for i, line in enumerate(read_PDB_file):
                    if 'CRYST1' in line:
                        read_x_dim_box_1 = read_PDB_file[i].split()[1:2]
                        read_x_dim_box_1 = float(read_x_dim_box_1[0])
                        read_x_dim_box_1_origin = read_x_dim_box_1 / 2
                        read_y_dim_box_1 = read_PDB_file[i].split()[2:3]
                        read_y_dim_box_1 = float(read_y_dim_box_1[0])
                        read_y_dim_box_1_origin = read_y_dim_box_1 / 2
                        read_z_dim_box_1 = read_PDB_file[i].split()[3:4]
                        read_z_dim_box_1 = float(read_z_dim_box_1[0])
                        read_z_dim_box_1_origin = read_z_dim_box_1 / 2

                        read_angle_alpha = read_PDB_file[i].split()[4:5]
                        read_angle_alpha = read_angle_alpha[0]
                        read_angle_beta = read_PDB_file[i].split()[5:6]
                        read_angle_beta = read_angle_beta[0]
                        read_angle_gamma = read_PDB_file[i].split()[6:7]
                        read_angle_gamma = read_angle_gamma[0]

                GOMC_used_x_dim = check_for_PDB_dims_and_override("x", Run_No, read_x_dim_box_1,
                                                                  set_dim = set_x_dim_box_1,
                                                                  only_on_Run_No = 1)
                GOMC_used_y_dim = check_for_PDB_dims_and_override("y", Run_No, read_y_dim_box_1,
                                                                  set_dim = set_y_dim_box_1,
                                                                  only_on_Run_No = 1)
                GOMC_used_z_dim = check_for_PDB_dims_and_override("z", Run_No, read_z_dim_box_1,
                                                                  set_dim=set_z_dim_box_1,
                                                                  only_on_Run_No = 1)

                new_GOMC_data = new_GOMC_data.replace("x_dim_box_1", str(GOMC_used_x_dim))
                new_GOMC_data = new_GOMC_data.replace("y_dim_box_1", str(GOMC_used_y_dim))
                new_GOMC_data = new_GOMC_data.replace("z_dim_box_1", str(GOMC_used_z_dim))

            else:
                new_GOMC_data = new_GOMC_data.replace("coor_box_1_file",
                                                      previous_GOMC_dir + "/" + "Output_data_BOX_1_restart.coor")
                new_GOMC_data = new_GOMC_data.replace("xsc_box_1_file",
                                                      previous_GOMC_dir + "/" + "Output_data_BOX_1_restart.xsc")

                previous_GOMC_xsc_box_1_file = previous_GOMC_dir + "/Output_data_BOX_1_restart.xsc"
                read_previous_GOMC_xsc_box_1_file = open(previous_GOMC_xsc_box_1_file, 'r').readlines()

                read_x_dim_box_1 = read_previous_GOMC_xsc_box_1_file[-1].split()[1:2]
                read_x_dim_box_1 = float(read_x_dim_box_1[0])
                read_y_dim_box_1 = read_previous_GOMC_xsc_box_1_file[-1].split()[5:6]
                read_y_dim_box_1 = float(read_y_dim_box_1[0])
                read_z_dim_box_1 = read_previous_GOMC_xsc_box_1_file[-1].split()[9:10]
                read_z_dim_box_1 = float(read_z_dim_box_1[0])

                read_x_dim_origin_box_1 = read_previous_GOMC_xsc_box_1_file[-1].split()[10:11]
                read_x_dim_origin_box_1 = float(read_x_dim_origin_box_1[0])
                read_y_dim_origin_box_1 = read_previous_GOMC_xsc_box_1_file[-1].split()[11:12]
                read_y_dim_origin_box_1 = float(read_y_dim_origin_box_1[0])
                read_z_dim_origin_box_1 = read_previous_GOMC_xsc_box_1_file[-1].split()[12:13]
                read_z_dim_origin_box_1 = float(read_z_dim_origin_box_1[0])

                new_GOMC_data = new_GOMC_data.replace("x_dim_box_1", str(read_x_dim_box_1))
                new_GOMC_data = new_GOMC_data.replace("y_dim_box_1", str(read_y_dim_box_1))
                new_GOMC_data = new_GOMC_data.replace("z_dim_box_1", str(read_z_dim_box_1))

        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
            NAMD_xsc_box_1_file = NAMD_box_1_newdir + "/namdOut.xsc"
            read_NAMD_xsc_box_1_file = open(NAMD_xsc_box_1_file, 'r').readlines()
            new_GOMC_data = new_GOMC_data.replace("coor_box_1_file", NAMD_box_1_newdir + '/namdOut.coor')
            new_GOMC_data = new_GOMC_data.replace("xsc_box_1_file", NAMD_box_1_newdir + '/namdOut.xsc')

            read_x_dim_box_1 = read_NAMD_xsc_box_1_file[-1].split()[1:2]
            read_x_dim_box_1 = float(read_x_dim_box_1[0])
            read_y_dim_box_1 = read_NAMD_xsc_box_1_file[-1].split()[5:6]
            read_y_dim_box_1 = float(read_y_dim_box_1[0])
            read_z_dim_box_1 = read_NAMD_xsc_box_1_file[-1].split()[9:10]
            read_z_dim_box_1 = float(read_z_dim_box_1[0])

            read_x_dim_origin_box_1 = read_NAMD_xsc_box_1_file[-1].split()[10:11]
            read_x_dim_origin_box_1 = float(read_x_dim_origin_box_1[0])
            read_y_dim_origin_box_1 = read_NAMD_xsc_box_1_file[-1].split()[11:12]
            read_y_dim_origin_box_1 = float(read_y_dim_origin_box_1[0])
            read_z_dim_origin_box_1 = read_NAMD_xsc_box_1_file[-1].split()[12:13]
            read_z_dim_origin_box_1 = float(read_z_dim_origin_box_1[0])

            new_GOMC_data = new_GOMC_data.replace("x_dim_box_1", str(read_x_dim_box_1))
            new_GOMC_data = new_GOMC_data.replace("y_dim_box_1", str(read_y_dim_box_1))
            new_GOMC_data = new_GOMC_data.replace("z_dim_box_1", str(read_z_dim_box_1))



    if simulation_type in ["GEMC", "GCMC"]:
        if simulation_type in ["GCMC"] and previous_GOMC_dir == "NA":
            new_GOMC_data = new_GOMC_data.replace("restart_true_or_false", 'false')
        elif simulation_type in [
            "GEMC"] and only_use_box_0_for_NAMD_for_GEMC == True and previous_GOMC_dir == "NA":
            new_GOMC_data = new_GOMC_data.replace("restart_true_or_false", 'false')
        else:
            new_GOMC_data = new_GOMC_data.replace("restart_true_or_false", 'true')
    else:
        new_GOMC_data = new_GOMC_data.replace("restart_true_or_false", 'true')

    new_GOMC_data = new_GOMC_data.replace("GOMC_Run_Steps", str((int(GOMC_Run_Steps))))
    new_GOMC_data = new_GOMC_data.replace("GOMC_RST_Coor_CKpoint_Steps",
                                          str((int(GOMC_RST_Coor_CKpoint_Steps))))
    new_GOMC_data = new_GOMC_data.replace("GOMC_console_BLKavg_Hist_Steps",
                                          str((int(GOMC_console_BLKavg_Hist_Steps))))
    new_GOMC_data = new_GOMC_data.replace("GOMC_Hist_sample_Steps", str(GOMC_Hist_sample_Steps))
    new_GOMC_data = new_GOMC_data.replace("System_temp_set", str(simulation_Temp_K))
    new_GOMC_data = new_GOMC_data.replace("System_press_set", str(simulation_Pressure_bar))

    if simulation_type in ["GCMC"]:
        if GCMC_only_ChemPot_Type_1_Residue == None and GCMC_only_ChemPot_Type_1_mu_K == None:
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_1", \
                                                  str(''))
        elif GCMC_only_ChemPot_Type_1_Residue != None and GCMC_only_ChemPot_Type_1_mu_K != None:
            print('GCMC_only_ChemPot_Type_1_mu_K = ' + str(GCMC_only_ChemPot_Type_1_mu_K))
            ChemPot_Type_1_string = str("ChemPot" + " \t " + GCMC_only_ChemPot_Type_1_Residue) \
                                    + " \t " + str(GCMC_only_ChemPot_Type_1_mu_K)
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_1", str(ChemPot_Type_1_string))
        else:
            write_log_data = "Warning: There is in error in the chemical potential settings for GCMC simulaiton" \
                             + ' \n'
            Log_Template_file.write(str(write_log_data))
            print(str(write_log_data))
            warn(str(write_log_data))

        if GCMC_only_ChemPot_Type_2_Residue == None and GCMC_only_ChemPot_Type_2_mu_K == None:
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_2", \
                                                  str(''))
        elif GCMC_only_ChemPot_Type_2_Residue != None and GCMC_only_ChemPot_Type_2_mu_K != None:
            ChemPot_Type_2_string = str("ChemPot" + " \t " + GCMC_only_ChemPot_Type_2_Residue) \
                                    + " \t " + str(GCMC_only_ChemPot_Type_2_mu_K)
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_2", str(ChemPot_Type_2_string))
        else:
            write_log_data = "Warning: There is in error in the chemical potential settings for GCMC simulaiton" \
                             + ' \n'
            Log_Template_file.write(str(write_log_data))
            print(str(write_log_data))
            warn(str(write_log_data))

        if GCMC_only_ChemPot_Type_3_Residue == None and GCMC_only_ChemPot_Type_3_mu_K == None:
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_3", \
                                                  str(''))
        elif GCMC_only_ChemPot_Type_3_Residue != None and GCMC_only_ChemPot_Type_3_mu_K != None:
            ChemPot_Type_3_string = str("ChemPot" + " \t " + GCMC_only_ChemPot_Type_3_Residue) \
                                    + " \t " + str(GCMC_only_ChemPot_Type_3_mu_K)
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_3", str(ChemPot_Type_3_string))
        else:
            write_log_data = "Warning: There is in error in the chemical potential settings for GCMC simulaiton" \
                             + ' \n'
            Log_Template_file.write(str(write_log_data))
            print(str(write_log_data))
            warn(str(write_log_data))

        if GCMC_only_ChemPot_Type_4_Residue == None and GCMC_only_ChemPot_Type_4_mu_K == None:
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_4", \
                                                  str(''))
        elif GCMC_only_ChemPot_Type_4_Residue != None and GCMC_only_ChemPot_Type_4_mu_K != None:
            ChemPot_Type_4_string = str("ChemPot" + " \t " + GCMC_only_ChemPot_Type_4_Residue) \
                                    + " \t " + str(GCMC_only_ChemPot_Type_4_mu_K)
            new_GOMC_data = new_GOMC_data.replace("mu_chempot_4", str(ChemPot_Type_4_string))
        else:
            write_log_data = "Warning: There is in error in the chemical potential settings for GCMC simulaiton" \
                             + ' \n'
            Log_Template_file.write(str(write_log_data))
            print(str(write_log_data))
            warn(str(write_log_data))

    try:
        new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(GOMC_Adj_Steps)))

    except:
        set_max_steps_equib_adj = 10 * 10 ** 6
        if GOMC_Run_Steps >= set_max_steps_equib_adj:
            new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int((set_max_steps_equib_adj / 10))))

            if GOMC_Run_Steps / 10 >= 1000:
                new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(1000)))
            else:
                new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(GOMC_Run_Steps / 10)))

        elif int(GOMC_Run_Steps / 10) > 0:
            # make equal to 1000 until the restart true can be enabled
            # new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(GOMC_Run_Steps/10)))
            new_GOMC_data = new_GOMC_data.replace("GOMC_Equilb_Steps", str((int(GOMC_Run_Steps / 10))))
            new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(GOMC_Run_Steps / 10)))

            if GOMC_Run_Steps / 10 >= 1000:
                new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(1000)))
            else:
                new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(GOMC_Run_Steps / 10)))

        else:
            # make equal to 1000 until the restart true can be enabled
            new_GOMC_data = new_GOMC_data.replace("GOMC_Adj_Steps", str(int(1)))
            new_GOMC_data = new_GOMC_data.replace("GOMC_Equilb_Steps", str((int(1))))

    if simulation_type in ["GEMC", "GCMC"]:
        if previous_GOMC_dir == 'NA':
            # marked as "Restart_Checkpoint_file", 'false' for now until checkpoint is setup
            new_GOMC_data = new_GOMC_data.replace("Restart_Checkpoint_file", 'false')
            new_GOMC_data = new_GOMC_data.replace("pdb_file_box_1_file",
                                                  str(Python_file_directory) + "/" + Starting_PDB_box_1_file)
            new_GOMC_data = new_GOMC_data.replace("psf_file_box_1_file",
                                                  str(Python_file_directory) + "/" + Starting_PSF_box_1_file)
        else:
            new_GOMC_data = new_GOMC_data.replace("pdb_file_box_1_file",
                                                  str(
                                                      previous_GOMC_dir) + "/" + "Output_data_BOX_1_restart.pdb")
            new_GOMC_data = new_GOMC_data.replace("psf_file_box_1_file",
                                                  str(
                                                      previous_GOMC_dir) + "/" + "Output_data_BOX_1_restart.psf")

    else:
        if previous_GOMC_dir == 'NA':
            # marked as "Restart_Checkpoint_file", 'false' for now until checkpoint is setup
            new_GOMC_data = new_GOMC_data.replace("Restart_Checkpoint_file", 'false')

    # make checkpoint true and restart
    new_GOMC_data = new_GOMC_data.replace("Restart_Checkpoint_file", str(use_GOMC_checkpoint))

    run_GOMC_copy_CkPt_new_dir_command = "ln -s " + str(previous_GOMC_dir) + "/" + "checkpoint.dat" + " " \
                                         + str(GOMC_newdir)

    exec_GOMC_copy_CkPt_new_dir_command = subprocess.Popen(run_GOMC_copy_CkPt_new_dir_command,
                                                           shell=True, stderr=subprocess.STDOUT)

    GOMC_copy_CkPt_pid_status = os.waitpid(exec_GOMC_copy_CkPt_new_dir_command.pid,
                                           os.WSTOPPED)  # pauses python until GOMC sim done

    generate_GOMC_file.write(new_GOMC_data)
    generate_GOMC_file.close()
    # *************************************************
    # build input file from template the GOMC simulation (end)
    # *************************************************

    return GOMC_newdir


def get_GOMC_energy_data(read_GOMC_box_X_log_file, box_number):
    # note GOMC energy units in kcal/mol
    # generate energy and system file data for box X
    Get_E_titles = True
    Get_V_P_TotMol_Rho_titles = True
    E_values_GOMC_List = []
    V_P_TotMol_Rho_values_List = []
    for i, line in enumerate(read_GOMC_box_X_log_file):
        if line.startswith('ETITLE:') == True and Get_E_titles == True:
            E_titles_GOMC_iteration = read_GOMC_box_X_log_file[i]
            E_titles_GOMC_iteration = E_titles_GOMC_iteration.split()
            E_titles_GOMC_iteration = E_titles_GOMC_iteration[:]
            if Get_E_titles == True:
                Get_E_titles = False

        if line.startswith('ENER_' + str(box_number)+ ':') == True:
            E_values_GOMC_iteration_List = []
            E_values_GOMC_iteration = read_GOMC_box_X_log_file[i]
            E_values_GOMC_iteration = E_values_GOMC_iteration.split()

            for j in range(0, len(E_values_GOMC_iteration)):
                if E_titles_GOMC_iteration[j] == 'ETITLE:':
                    E_values_GOMC_iteration_List.append(E_values_GOMC_iteration[j])
                elif E_titles_GOMC_iteration[j] == 'STEP':
                    E_values_GOMC_iteration_List.append(int(int(E_values_GOMC_iteration[j]) + current_step))
                else:
                    E_values_GOMC_iteration_List.append(np.float(E_values_GOMC_iteration[j]) * K_to_kcal_mol)
            E_values_GOMC_List.append(E_values_GOMC_iteration_List)

        if line.startswith('STITLE:') == True and Get_V_P_TotMol_Rho_titles == True:
            V_P_TotMol_Rho_titles_iteration = read_GOMC_box_X_log_file[i]
            V_P_TotMol_Rho_titles_iteration = V_P_TotMol_Rho_titles_iteration.split()
            V_P_TotMol_Rho_titles_iteration = V_P_TotMol_Rho_titles_iteration[:]
            if Get_V_P_TotMol_Rho_titles == True:
                Get_V_P_TotMol_Rho_titles = False

        if line.startswith('STAT_' + str(box_number)+ ':') == True:
            V_P_TotMol_Rho_values_iteration_List = []
            V_P_TotMol_Rho_values_iteration = read_GOMC_box_X_log_file[i]
            V_P_TotMol_Rho_values_iteration = V_P_TotMol_Rho_values_iteration.split()

            for j in range(0, len(V_P_TotMol_Rho_values_iteration)):
                if V_P_TotMol_Rho_titles_iteration[j] == 'STITLE:':
                    V_P_TotMol_Rho_values_iteration_List.append(V_P_TotMol_Rho_values_iteration[j])
                elif V_P_TotMol_Rho_titles_iteration[j] == 'STEP':
                    V_P_TotMol_Rho_values_iteration_List.append(int(int(V_P_TotMol_Rho_values_iteration[j]) \
                                                                    + current_step))
                else:
                    V_P_TotMol_Rho_values_iteration_List.append(np.float(V_P_TotMol_Rho_values_iteration[j]))
            V_P_TotMol_Rho_values_List.append(V_P_TotMol_Rho_values_iteration_List)

    GOMC_energy_data_box_X_df = pd.DataFrame(data=E_values_GOMC_List, columns=E_titles_GOMC_iteration)

    return GOMC_energy_data_box_X_df

def get_GOMC_energy_data_kcal_per_mol(GOMC_energy_data_box_X_df):
    GOMC_E_electro_box_X_kcal_per_mol = GOMC_energy_data_box_X_df.loc[:, 'TOTAL_ELECT'].tolist()
    GOMC_E_electro_box_X_initial_value = GOMC_E_electro_box_X_kcal_per_mol[0]
    GOMC_E_electro_box_X_final_value = GOMC_E_electro_box_X_kcal_per_mol[-1]

    GOMC_E_potential_box_X_kcal_per_mol = GOMC_energy_data_box_X_df.loc[:, 'TOTAL'].tolist()
    GOMC_E_potential_box_X_initial_value = GOMC_E_potential_box_X_kcal_per_mol[0]
    GOMC_E_potential_box_X_final_value = GOMC_E_potential_box_X_kcal_per_mol[-1]

    GOMC_E_INTRA_NB_box_X_kcal_per_mol = GOMC_energy_data_box_X_df.loc[:, 'INTRA(NB)'].tolist()
    GOMC_E_INTRA_NB_box_X_initial_value = GOMC_E_INTRA_NB_box_X_kcal_per_mol[0]
    GOMC_E_INTRA_NB_box_X_final_value = GOMC_E_INTRA_NB_box_X_kcal_per_mol[-1]

    GOMC_E_INTER_LJ_box_X_kcal_per_mol = GOMC_energy_data_box_X_df.loc[:, 'INTER(LJ)'].tolist()
    GOMC_E_INTER_LJ_box_X_initial_value = GOMC_E_INTER_LJ_box_X_kcal_per_mol[0]
    GOMC_E_INTER_LJ_box_X_final_value = GOMC_E_INTER_LJ_box_X_kcal_per_mol[-1]

    GOMC_E_vdw_box_X_kcal_per_mol = []
    for vwd_i in range(0, len(GOMC_E_INTRA_NB_box_X_kcal_per_mol)):
        GOMC_E_vdw_box_X_kcal_per_mol.append(GOMC_E_INTRA_NB_box_X_kcal_per_mol[vwd_i] + \
                                             GOMC_E_INTER_LJ_box_X_kcal_per_mol[vwd_i])
    GOMC_E_vdw_box_X_initial_value = GOMC_E_vdw_box_X_kcal_per_mol[0]
    GOMC_E_vdw_box_X_final_value = GOMC_E_vdw_box_X_kcal_per_mol[-1]

    return GOMC_E_electro_box_X_kcal_per_mol, GOMC_E_electro_box_X_initial_value, \
           GOMC_E_electro_box_X_final_value, \
           GOMC_E_potential_box_X_kcal_per_mol, GOMC_E_potential_box_X_initial_value, \
           GOMC_E_potential_box_X_final_value, \
           GOMC_E_vdw_box_X_kcal_per_mol, GOMC_E_vdw_box_X_initial_value, GOMC_E_vdw_box_X_final_value







for Run_No in range(Starting_sims_NAMD_GOMC, Total_sims_NAMD_GOMC):
    # *************************************************
    # *************************************************
    # Simulation initial or restart setup (Start)
    # *************************************************
    # *************************************************
    #set_box_numbers
    box_number_0 = 0
    box_number_1 = 1

    try:
        print("current step = " +str(int(current_step- NAMD_Minimize_Steps)) + " ( negative (-) for minimize steps)")
    except:
        print("current step = NA " + " ( negative (-) for minimize steps)")

    add_zeros_at_start_Run_No_str = Calc_folder_zeros(Run_No)

    # setting the past file name in a restart
    if Run_No == Starting_sims_NAMD_GOMC and Starting_sims_NAMD_GOMC != 0:
        add_zeros_at_start_Run_No_str = Calc_folder_zeros(int(Starting_sims_NAMD_GOMC - 2))
        NAMD_box_0_newdir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                            str(add_zeros_at_start_Run_No_str) + str(int(Starting_sims_NAMD_GOMC - 2)) + "_a"

        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
            NAMD_box_1_newdir = str(Python_file_directory) + "/" + path_NAMD_runs + "/" + \
                                str(add_zeros_at_start_Run_No_str) + str(int(Starting_sims_NAMD_GOMC - 2)) + "_b"

        add_zeros_at_start_Run_No_str = Calc_folder_zeros(int(Starting_sims_NAMD_GOMC -1))
        GOMC_newdir = str(Python_file_directory) + "/" + path_GOMC_runs + "/" \
                      + str(add_zeros_at_start_Run_No_str) + str(int(Starting_sims_NAMD_GOMC -1))

        current_step = (NAMD_Run_Steps + GOMC_Run_Steps) * Starting_at_cycle_NAMD_GOMC_sims + NAMD_Minimize_Steps  # steps in number of cycles

    elif Run_No == Starting_sims_NAMD_GOMC and Starting_sims_NAMD_GOMC == 0:
        current_step = 0

    # Delete FFT info from NAMD if starting a new simulation (i.e., Run_No at start of running == 0)
    if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False and Run_No == 0:
        delete_NAMD_run_0_FFT_file(box_number_0)
        delete_NAMD_run_0_FFT_file(box_number_1)

    elif Run_No == 0:
        delete_NAMD_run_0_FFT_file(box_number_0)

    # Get the NAMD_PME_GRID_DIMs upon restarting a simulation
    if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False and \
            (Run_No == Starting_sims_NAMD_GOMC ):
        NAMD_X_PME_GRID_box_0_DIM, NAMD_Y_PME_GRID_box_0_DIM, \
        NAMD_Z_PME_GRID_box_0_DIM, NAMD_box_0_Run_0_dir = get_NAMD_run_0_PME_DIM(box_number_0)
        NAMD_X_PME_GRID_box_1_DIM, NAMD_Y_PME_GRID_box_1_DIM, \
        NAMD_Z_PME_GRID_box_1_DIM, NAMD_box_1_Run_0_dir = get_NAMD_run_0_PME_DIM(box_number_1)

    elif Starting_sims_NAMD_GOMC == Run_No and  (Run_No == Starting_sims_NAMD_GOMC):
        NAMD_X_PME_GRID_box_0_DIM, NAMD_Y_PME_GRID_box_0_DIM, \
        NAMD_Z_PME_GRID_box_0_DIM, NAMD_box_0_Run_0_dir = get_NAMD_run_0_PME_DIM(box_number_0)

    write_log_data = "*************************************************\n" + \
                     "*************************************************\n" + \
                     'Run_No = '+ str(Run_No ) + " (START)\n"\
                     "*************************************************" + ' \n'
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))

    if Run_No != 0  :
        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False \
                and Starting_sims_NAMD_GOMC == Run_No:
            NAMD_box_0_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number_0)
            NAMD_box_1_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number_1)

        else:
            NAMD_box_0_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number_0)
    # *************************************************
    # *************************************************
    # Simulation initial or restart setup (end)
    # *************************************************
    # *************************************************


    # *************************************************
    # *************************************************
    # RUN THE NAMD PORTION of the CODE (Start)
    # *************************************************
    # *************************************************
    if Run_No % 2 == 0:  # NAMD's run time. :# NAMD's staring run.  NAMD starts simulation series with energy minimization
        # *************************************************
        # build input file from template for box 0 (start)
        # *************************************************

        try:
            GOMC_newdir
        except:
            GOMC_newdir = "NA"

        NAMD_box_0_newdir = write_NAMD_conf_file(Python_file_directory, path_NAMD_template,  path_NAMD_runs,
                                                 GOMC_newdir, Run_No, box_number_0, NAMD_Run_Steps,
                                                 NAMD_Minimize_Steps, NAMD_RST_DCD_XST_Steps,
                                                 NAMD_console_BLKavg_E_and_P_Steps,
                                                 simulation_Temp_K, simulation_Pressure_bar,
                                                 Starting_PDB_box_0_file, Starting_PSF_box_0_file,
                                                 NAMD_X_PME_GRID_box_0_DIM,
                                                 NAMD_Y_PME_GRID_box_0_DIM,
                                                 NAMD_Z_PME_GRID_box_0_DIM,
                                                 set_x_dim = set_x_dim_box_0,
                                                 set_y_dim = set_y_dim_box_0,
                                                 set_z_dim = set_z_dim_box_0)

        # *************************************************
        # build input file from template for box 0 (end)
        # *************************************************

        # *************************************************
        # build input file from template for box 1 (start)
        # *************************************************
        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:

            NAMD_box_1_newdir = write_NAMD_conf_file(Python_file_directory, path_NAMD_template, path_NAMD_runs,
                                                     GOMC_newdir, Run_No, box_number_1, NAMD_Run_Steps,
                                                     NAMD_Minimize_Steps, NAMD_RST_DCD_XST_Steps,
                                                     NAMD_console_BLKavg_E_and_P_Steps,
                                                     simulation_Temp_K, simulation_Pressure_bar,
                                                     Starting_PDB_box_1_file, Starting_PSF_box_1_file,
                                                     NAMD_X_PME_GRID_box_1_DIM,
                                                     NAMD_Y_PME_GRID_box_1_DIM,
                                                     NAMD_Z_PME_GRID_box_1_DIM,
                                                     set_x_dim = set_x_dim_box_1,
                                                     set_y_dim = set_y_dim_box_1,
                                                     set_z_dim = set_z_dim_box_1,
                                                     FFT_add_NAMD_Ang_to_box_dim=0)
            # *************************************************
            # build input file from template for box 1 (end)
            # *************************************************

        # *************************************************
        # Run the intitial NAMD simulations box 0 and 1 (start)
        # *************************************************
        write_log_data = "*************************************************\n" \
                         +'Runnging the initial NAMD simulations now' + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

        if Run_No != 0:
            # copy the FFT grid.txt file from the first NAMD simulation (i.e., Run 0) to the current dir
            # for box_0
            cp_NAMD_box_0_FFT_Run_0_new_dir_cmd = "ln -s " + str(NAMD_box_0_Run_0_dir) + "/" \
                                                  + str(NAMD_box_0_Run_0_FFT_NAMD_filename) + " " \
                                                  + str(NAMD_box_0_newdir)

            exec_NAMD_box_0_cp_FFT_Run_0_new_dir_cmd = subprocess.Popen(cp_NAMD_box_0_FFT_Run_0_new_dir_cmd,
                                                                        shell=True, stderr=subprocess.STDOUT)

        run_box_0_command = "cd " + str(NAMD_box_0_newdir) + " && " \
                            + str(NAMD_executable_file) + " +p" \
                            + str(int(Total_No_cores)) + " in.conf > out.dat" + " "

        exec_run_box_0_command = subprocess.Popen(run_box_0_command, shell=True, stderr=subprocess.STDOUT)

        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
            # find the FFT grid.txt file from the first NAMD simulation (i.e., Run 0) and copy it to
            # to the current dir for box 0
            if Run_No != 0:
                # copy the FFT grid.txt file from the first NAMD simulation (i.e., Run 0) to the current dir
                # for box_0 and box_1

                cp_NAMD_box_1_FFT_Run_0_new_dir_cmd = "cp " + str( NAMD_box_1_Run_0_dir) + "/" \
                                                             + str(NAMD_box_1_Run_0_FFT_NAMD_filename) + " " \
                                                             + str(NAMD_box_1_newdir)

                exec_NAMD_box1_cp_FFT_Run_0_new_dir_cmd = subprocess.Popen(cp_NAMD_box_1_FFT_Run_0_new_dir_cmd,
                                                                           shell=True, stderr=subprocess.STDOUT)

            run_box_1_command = "cd " + str(NAMD_box_1_newdir) + " && " \
                                + str(NAMD_executable_file) + " +p" \
                                + str(int(No_core_box_1)) + " in.conf > out.dat" + " "

            exec_run_box_1_command = subprocess.Popen(run_box_1_command, shell=True, stderr=subprocess.STDOUT)

        write_log_data = 'Waiting for initial NAMD simulation to finish ' + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

        box_0_pid_status = os.waitpid(exec_run_box_0_command.pid, os.WSTOPPED)  # pauses python until box 0 sim done
        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
            box_1_pid_status = os.waitpid(exec_run_box_1_command.pid, os.WSTOPPED) # pauses python until box 0 sim done

        write_log_data = 'The NAMD simulation are finished ' + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))


        # *************************************************
        # Run the intitial NAMD simulations box 0 and 1 (end)
        # *************************************************


        #*************************************
        # get final system energies for box 0 and 1 (start)
        # ***********************initial_Energies**************
        read_NAMD_box_0_Energy_file = open(str(NAMD_box_0_newdir) + "/" + 'out.dat', 'r').readlines()
        # note NAMD energy units in kcal/mol (no modifications required)
        # generate energy file data for box 0

        NAMD_E_electro_box_0, \
        NAMD_E_electro_box_0_initial_value, \
        NAMD_E_electro_box_0_final_value, \
        NAMD_E_potential_box_0, \
        NAMD_E_potential_box_0_initial_value, \
        NAMD_E_potential_box_0_final_value, \
        NAMD_E_vdw_box_0, \
        NAMD_E_vdw_box_0_initial_value, \
        NAMD_E_vdw_box_0_final_value = get_NAMD_energy_data(read_NAMD_box_0_Energy_file,
                                                            default_NAMD_E_titles )


        # note NAMD energy units in kcal/mol (no modifications required)
        # generate energy file data for box 1
        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
            read_NAMD_box_1_Energy_file = open(str(NAMD_box_1_newdir) + "/" + 'out.dat', 'r').readlines()

            NAMD_E_electro_box_1, \
            NAMD_E_electro_box_1_initial_value, \
            NAMD_E_electro_box_1_final_value, \
            NAMD_E_potential_box_1, \
            NAMD_E_potential_box_1_initial_value, \
            NAMD_E_potential_box_1_final_value, \
            NAMD_E_vdw_box_1, \
            NAMD_E_vdw_box_1_initial_value, \
            NAMD_E_vdw_box_1_final_value = get_NAMD_energy_data(read_NAMD_box_1_Energy_file,
                                                                default_NAMD_E_titles)


        if Run_No != 0 and Run_No != Starting_sims_NAMD_GOMC:
            # Compare the Last GOMC and first NAMD value to confirm the simulation data
            # VMD comparison between NAMD and GOMC data box 0

            compare_NAMD_GOMC_energies(GOMC_E_vdw_box_0_final_value, NAMD_E_vdw_box_0_initial_value,
                                       GOMC_E_electro_box_0_final_value, NAMD_E_electro_box_0_initial_value,
                                       Run_No, box_number_0)
            if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
                compare_NAMD_GOMC_energies(GOMC_E_vdw_box_1_final_value, NAMD_E_vdw_box_1_initial_value,
                                           GOMC_E_electro_box_1_final_value, NAMD_E_electro_box_1_initial_value,
                                           Run_No, box_number_1)

        # get the NAMD FFT file name to copy for future NAMD simulations if starting a new NAMD/GOMC simulation
        if Run_No == 0:
            if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False \
                    and Starting_sims_NAMD_GOMC == Run_No:
                NAMD_box_0_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number_0)
                NAMD_box_1_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number_1)
            else:
                NAMD_box_0_Run_0_FFT_NAMD_filename = get_NAMD_run_0_FFT_filename(box_number_0)

        #Get the NAMD_PME_GRID_DIMs upon finishing a new simulation
        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False and \
                (Run_No == Starting_sims_NAMD_GOMC ):
            NAMD_X_PME_GRID_box_0_DIM, NAMD_Y_PME_GRID_box_0_DIM, \
            NAMD_Z_PME_GRID_box_0_DIM, NAMD_box_0_Run_0_dir = get_NAMD_run_0_PME_DIM(box_number_0)
            NAMD_X_PME_GRID_box_1_DIM, NAMD_Y_PME_GRID_box_1_DIM, \
            NAMD_Z_PME_GRID_box_1_DIM, NAMD_box_1_Run_0_dir = get_NAMD_run_0_PME_DIM(box_number_1)

        elif Starting_sims_NAMD_GOMC == Run_No and (Run_No == Starting_sims_NAMD_GOMC):
            NAMD_X_PME_GRID_box_0_DIM, NAMD_Y_PME_GRID_box_0_DIM, \
            NAMD_Z_PME_GRID_box_0_DIM, NAMD_box_0_Run_0_dir = get_NAMD_run_0_PME_DIM(box_number_0)


        if Run_No != 0:
            current_step += NAMD_Run_Steps
        else:
            current_step += NAMD_Run_Steps + NAMD_Minimize_Steps
        # *************************************
        # get final system energies for box 0 and 1 (start)
        # *************************************
    # *************************************************
    # *************************************************
    # RUN THE NAMD PORTION of the CODE (End)
    # *************************************************
    # *************************************************






    # *************************************************
    # *************************************************
    # RUN THE GOMC PORTION of the CODE (Start)
    # *************************************************
    # *************************************************
    elif Run_No % 2 == 1:  # GOMC's run time.  GOMC starts simulation series
        # GOMC runs
        # *************************************************
        # build input file from template the GOMC simulation (start)
        # *************************************************
        try:
            previous_GOMC_dir = GOMC_newdir
        except:
            previous_GOMC_dir = "NA"

        GOMC_newdir =  write_GOMC_conf_file(Python_file_directory, path_GOMC_runs, Run_No, GOMC_Run_Steps,
                                            GOMC_RST_Coor_CKpoint_Steps, GOMC_console_BLKavg_Hist_Steps,
                                            GOMC_Hist_sample_Steps,
                                            simulation_Temp_K, simulation_Pressure_bar,
                                            Starting_PDB_box_0_file, Starting_PDB_box_1_file)

        write_log_data = "GOMC simulation data for simulation number " + str(Run_No) + " is completed" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))
        # *************************************************
        # build input file from template the GOMC simulation (start)
        # *************************************************

        # *************************************************
        # Copy file and Run the GOMC simulations (start)
        # *************************************************

        write_log_data = "*************************************************\n" \
                         + 'Runnging the GOMC simulations now ' + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))
        run_GOMC_command = "cd " + str(GOMC_newdir) + " && " \
                            + str(GOMC_executable_file) + " +p" \
                            + str(int(Total_No_cores)) + " in.conf > out.dat" + " "


        exec_GOMC_run_command = subprocess.Popen(run_GOMC_command, shell=True, stderr=subprocess.STDOUT)

        write_log_data = 'Waiting for initial GOMC simulation to finish '
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))
        GOMC_pid_status = os.waitpid(exec_GOMC_run_command.pid, os.WSTOPPED)  # pauses python until box 0 sim done
        write_log_data = 'The GOMC simulation are finished '
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))



        write_log_data = "Completed simulation in GOMC command\n" \
                         + "*************************************************" + ' \n'
        Log_Template_file.write(str(write_log_data))
        print(str(write_log_data))

        # *************************************************
        # Copy file and Run the GOMC simulations (end)
        # *************************************************



        #*******************************************************
        # get final system energies for box 0 and 1 (start)
        # ***********************initial_Energies**************
        read_GOMC_log_file = open(str(GOMC_newdir) + "/" + 'out.dat', 'r').readlines()
        GOMC_energy_data_box_0_df = get_GOMC_energy_data(read_GOMC_log_file, box_number_0)

        if simulation_type in ["GEMC", "GCMC"]:
            GOMC_energy_data_box_1_df = get_GOMC_energy_data(read_GOMC_log_file, box_number_1)


        # retrieve energy data from the printed file for the first and last points for box 0
        # extrated in units of kcal per mol
        GOMC_E_electro_box_0_kcal_per_mol, \
        GOMC_E_electro_box_0_initial_value, \
        GOMC_E_electro_box_0_final_value, \
        GOMC_E_potential_box_0_kcal_per_mol, \
        GOMC_E_potential_box_0_initial_value, \
        GOMC_E_potential_box_0_final_value, \
        GOMC_E_vdw_box_0_kcal_per_mol, \
        GOMC_E_vdw_box_0_initial_value, \
        GOMC_E_vdw_box_0_final_value = get_GOMC_energy_data_kcal_per_mol(GOMC_energy_data_box_0_df)

        # retrieve energy data from the printed file for the first and last points for box 1s
        if simulation_type in ["GEMC", "GCMC"]:  
            # extrated in units of kcal per mol
            GOMC_E_electro_box_1_kcal_per_mol, \
            GOMC_E_electro_box_1_initial_value, \
            GOMC_E_electro_box_1_final_value, \
            GOMC_E_potential_box_1_kcal_per_mol, \
            GOMC_E_potential_box_1_initial_value, \
            GOMC_E_potential_box_1_final_value, \
            GOMC_E_vdw_box_1_kcal_per_mol, \
            GOMC_E_vdw_box_1_initial_value, \
            GOMC_E_vdw_box_1_final_value = get_GOMC_energy_data_kcal_per_mol(GOMC_energy_data_box_1_df)


        # Compare the Last NAMD and first GOMC value to confirm the simulation data

        # VMD comparison between NAMD and GOMC data box 0
        compare_NAMD_GOMC_energies(NAMD_E_vdw_box_0_final_value, GOMC_E_vdw_box_0_initial_value,
                                   NAMD_E_electro_box_0_final_value, GOMC_E_electro_box_0_initial_value,
                                   Run_No, box_number_0)
        if simulation_type in ["GEMC"] and only_use_box_0_for_NAMD_for_GEMC == False:
            compare_NAMD_GOMC_energies(NAMD_E_vdw_box_1_final_value, GOMC_E_vdw_box_1_initial_value,
                                       NAMD_E_electro_box_1_final_value, GOMC_E_electro_box_1_initial_value,
                                       Run_No, box_number_1)


        current_step += GOMC_Run_Steps
        # *************************************************
        # build input file from template GOMC box 0 (end)
        # *************************************************

    # *************************************************
    # *************************************************
    # RUN THE NAMD PORTION of the CODE (End)
    # *************************************************
    # *************************************************


    write_log_data = "*************************************************\n" + \
                     'Run_No = '+ str(Run_No ) + " (End) "  + ' \n'
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))

end_time = datetime.datetime.today()
write_log_data = "*************************************************\n"  \
                 + 'date and time (end) = ' + str(end_time ) + " \n" \
                 + 'total simulation time = ' + str(end_time- start_time ) + " \n" \
                 + "*************************************************" + ' \n'
Log_Template_file.write(str(write_log_data))
print(str(write_log_data))

Log_Template_file.close()
