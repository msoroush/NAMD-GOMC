#*********************************
# TO RUN hybrid simulations: (Start)
#*********************************
1) Change the variables thru line ~80 in the run "run_NAMD_GOMC_Rev_X.py" file.  Note the information about the variables is in the python code.  

2) If you have previous  "NAMD" and  "GOMC" folders in this directory deleted them, so you do not have mixed data in the simulation folders.  
 
3) Assuming that you have the packages installed or running anaconda env, if not install them.  Then run :

python run_NAMD_GOMC_Rev_14.py       # runs the hybrid simulations (for now if its in its current directory)

4) The simulation runs are sent to the "NAMD" and  "GOMC" folders, and the hybrid log file is in the same directory as the  combined data is now located in the "run_NAMD_GOMC_Rev_X.py" file .
#*********************************
# TO RUN hybrid simulations: (end)
#*********************************



#*********************************
# Run the code to combine the data (i.e., combine_data_NAMD_GOMC_Rev_X.py code) (Start)
#*********************************
1) This file currently needs to be run from the directory with the "NAMD" and  "GOMC" folders !!!!

2) Change the variables to the same as the ones used in the code:
	-  simulation_type     			 (example:   simulation_type = 'GEMC' )
	-  only_use_box_0_for_NAMD_for_GEMC      (example:   only_use_box_0_for_NAMD_for_GEMC = False )
	-  Simulation_Engine_Options 		 (example:   Simulation_Engine_Options = 'Hybrid'  )
	-  GOMC_only_log_filename		 (example:   GOMC_only_log_filename = '../../MC_run_GEMC/out.dat' )

3) Assuming that you have the packages installed or running anaconda env, if not install them.  Then run :

python combine_data_NAMD_GOMC_Rev_6.py    # runs all the combining file data provided its in its current location

4) The combined data is now located in the "combined_data" folder, with the minimization run data removed (i.e., step 0 is the first step after the minimization.)
#*********************************
# Run the code to combine the data (i.e., combine_data_NAMD_GOMC_Rev_X.py code) (end)
#*********************************



#*********************************
# Plot the GOMC extracted data: total Potential energy and the density (Start)
#*********************************
1) This file currently the user to selected the directory with the "GOMC-only" and  "Hybrid"  combined data files (i.e., "GOMC_Energies_Stat_kcal_per_mol_box_0.txt'" files)!!!!

2) Change the variables thru line ~40 for your system. The Box_No variable is provided to use easy switching between box plotting.

3) Assuming that you have the packages installed or running anaconda env, if not install them.  Then run :

python Plot_box_Data_GOMC_only_vs_Hybrid_rev_1.py    # adds all the plots to the directory it is in


#*********************************
# Plot the GOMC extracted data: total Potential energy and the density (end)
#*********************************
