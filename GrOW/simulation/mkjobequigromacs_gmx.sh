#!/bin/bash

# ****************************************************************************
#
# mkjobequigromacs - Gromacs execution subroutine for GROW
# ========================================================
# 
# This subroutine runs a full set of simulations to get the results of a fully 
# equilibrated simulation.
#
# original code: M. Huelsmann
# modifications: T.J. Mueller (November 9th 2009)
#
# ****************************************************************************

# ****************************************************************************
#
# Structure of mkjobequigromacs.sh
# ================================
#
# 1. Preparations
#    1.1. Variables
#    1.2. File handling
#    1.3. Create new topology (Substitute place holders of input file)
#    1.4. Check output to see whether variables are transported correctly
#    1.5. Rerun option
# 2. Energy minimization
#    2.1. Main title
#    2.2. Make energy minimization directory
#    2.3. Preparation routines for simulation tool
#         2.3.1. Information
#         2.3.2. Execution
#    2.4. Run the simulation tool executable
#         2.4.1. Information
#         2.4.2. Execution
#    2.5. Return to the parent directory
# 3. Pre-pre-equilibration
#    3.1. Main title
#    3.2. Make pre-pre-equilibration directory
#    3.3. Update Parameters for the pre-pre-equilibration
#    3.4. Preparation routines for simulation tool
#         3.4.1. Information
#         3.4.2. Execution
#    3.5. Run the simulation tool executable
#         3.5.1. Information
#         3.5.2. Execution
#    3.6. Return to the parent directory
# 4. Pre-equilibration
#    4.1. Main title
#    4.2. Make pre-pre-equilibration directory
#    4.3. Update Parameters for the preequilibration
#         4.3.1. Set time step, number of time steps, and restart time, if necessary
#         4.3.2. Update ini file according to commandline options
#    4.4. Preparation routines for simulation tool
#         4.4.1. Information
#         4.4.2. Execution
#    4.5. Run the simulation tool executable
#         4.5.1. Information
#         4.5.2. Execution
#    4.6. Return to the parent directory
# 5. Equilibration cycles
#    5.1. Main title
#    5.2. Preparing the cycles
#         5.2.1 Create equilibration directory and go there
#         5.2.2 Redefine simulation variables
#         5.2.3 Create evaluation files for simulation tool, if necessary
#         5.2.4 data_all (Potential energy)
#    5.3. Start cycles
#         5.3.1 Modifying the simulation ini file
#         5.3.2 Title of cycle
#         5.3.3 Clean files
#         5.3.4 Preparation routines for simulation tool
#               5.3.4.1 Information
#               5.3.4.2 Execution
#         5.3.5 Run the simulation tool executable
#               5.3.5.1 Information
#               5.3.5.2 Execution
#         5.3.6 Evaluating the current cycle
#               5.3.6.1 Title
#               5.3.6.2 Get general property of interest (potential energy)
#               5.3.6.3 Get density
#               5.3.6.4 Get nonbonded energy
#               5.3.6.5 Convert files into GROW standard lists
#         5.3.7 Analyse the current data
#               5.3.7.1 Title
#               5.3.7.2 Append the evaluated data to a total equilibration evaluation file
#               5.3.7.3 Extract from all data the last $nof_samples data points
#               5.3.7.4 Do the statistics via R
#               5.3.7.5 Finishing equilibration: writing final coordinates file for production
#    5.4. Return to the parent directory
# 6. Production run
#    6.1. Main title
#    6.2. Make production directory
#    6.3. Modifying simulation ini file
#         6.3.1 Set time step, number of time steps, and restart time, if necessary
#         6.3.2 Update ini file according to commandline options
#    6.4. Preparation routines for simulation tool
#         6.4.1 Information
#         6.4.2 Execution
#    6.5. Run the simulation tool executable
#         6.5.1 Information
#         6.5.2 Execution
#    6.6. Evaluate production run
#         6.6.1 Calculate density, temperature, volume, pressure, and potential energy
#         6.6.2 Calculate diffusion coefficient
#         6.6.3 Calculate reorientation time
#         6.6.4 Calculate viscosity
#         6.6.5 Calculate nonbonded energy and enthalpy of vaporization
#         6.6.6 Calculate molecular weight
#         6.6.7 Calculate other properties required for optimization workflow
# 7. Clean up and Termination
#
# ****************************************************************************

# *****************************************************************************
# 1. Preparations
# *****************************************************************************

# 1.1 Variables
# =============

# Language
export LANG=en_US.UTF-8
export MPIRUN_TIMEOUT=36000
export PATH=$PATH:$BINDIR:$MOSCITOBINPATH:$FITPOLYBINPATH

# Verbosity
verbose=$VERBOSITY

# simulation specific variables
np=$NODES
gromacs_dir=$BINDIR
topology=$TPDUMMY
in_coordinates=$COORDS
sim_parameters=$equi_in
input_dir=$INPUTDIR
quantity_of_interest=Potential
r_script=$RSCRIPT
nof_steps=$NSTEPS
temperature=$TEMPERATURE
pressure=$PRESSURE
timestep=$TIMESTEP
restart_time=1000     

# equilibration specific variables
cycle=0

# 1.2 File handling
# =================

# set working path
working_dir=$OUTPATH
# alternatively: working_dir=`pwd`

# keep the current path in mind
act_dir=`pwd`

# initialize a new subdirectory
if ! test -e $working_dir/subdir; then
   mkdir $working_dir/subdir
fi

cd $working_dir/subdir

# 1.3. Create new topology (Substitute place holders of input file)
# =================================================================

# parameters in $@ have to be written into program specific input file containing
# place holders for the force field parameters to be optimized

head='sed '
substitute=''
tudels='" '
tail=' $topology > $working_dir/current.top'
index=1

for par in $@; do

	substitute="$substitute $sed s/X_$index/$par/;"
	index=`expr $index + 1` 


done

eval $head$tudels$substitute$tudels$tail

head='sed '
substitute=''
tudels='" '
tail=' $diff_system > $working_dir/current.sys'
index=1

echo $CPARSTRING

for par in $CPARSTRING; do

        substitute="$substitute $sed s/X_$index/$par/;"
        index=`expr $index + 1`


done

eval $head$tudels$substitute$tudels$tail



# 1.4. Check output to see whether variables are transported correctly
# ====================================================================

if [[ $verbose -gt 5 ]]; then
    echo "mkjobequigromacs.sh: check transfer of input"
    echo " "

    echo "     np                   .$np."
    echo "     gromacs_dir          .$gromacs_dir."	          
    echo "     input_dir            .$input_dir."   
    echo "     topology             .$topology."
    echo "     in_coordinates       .$in_coordinates."	     
    echo "     sim_parameters       .$sim_parameters."	     
    echo "     standard quantity .$quantity_of_interest."
    if [[ $DENSITY == "y" ]]; then echo "     quantity .Density."; fi
    if [[ $NONBOND == "y" ]]; then echo "     quantity .Nonbonded Energy."; fi
    if [[ $NONBOND == "y" ]]; then echo "vapor_script         .$VAPORSCRIPT."; fi
    echo "     r_script             .$r_script."
    echo "     NMOL                  .$NMOL."
    echo "     nof_steps            .$nof_steps."
    echo "     temperature          .$temperature."	     
    echo "     pressure             .$pressure."
    echo "     timestep             .$timestep."
    echo "     restart_time         .$restart_time."  
    echo "     output_dir           .$working_dir."	  
    echo " "
fi



# 1.5. Rerun option
# ====================================================================

rerun_exec=FALSE

if [[ $RERUN == "y" ]]; then
   string1=${working_dir/*\//}
   string=${string1/.out.*/}
   
   if test $string == gradient; then   
      rerun_exec=TRUE
   fi
fi

# =======================================================================

if [[ $rerun_exec == "FALSE" ]]; then

# *****************************************************************************
# 2. Energy minimization
# *****************************************************************************

# assuming the initial state being far away from equilibrium, the system may
# need a cautious start
# preparation and simulation must be provided for the simulation tool

# 2.1. Main title
# ===============

if [[ $verbose -gt 3 ]]; then
    echo " Energy minimization"
    echo " ==================="
    echo " "
fi

# 2.2. Make energy minimization directory
# =======================================

if ! test -e energy_minimization; then
    mkdir energy_minimization
fi

cd energy_minimization

# 2.3. Preparation routines for simulation tool
# =============================================

# must be exchanged for a different tool


# 2.3.1 Information
# -----------------

   if [[ $verbose -gt 5 ]]; then 
    echo " gmx grompp: Input md parameters (-f): " $EMIN
    echo "         Input topology (-p) :     " $working_dir/current.top
    echo "         input structure (-c):     " $in_coordinates
    echo "         standard output:          " e-min_gmx_grompp 
    echo "" 
   fi

# 2.3.2 Execution
# ---------------
   gmx grompp >& e-min_gmx_grompp \
        -f $EMIN \
        -p $working_dir/current.top \
        -o $working_dir/subdir/energy_minimization/e-min.tpr \
        -c $in_coordinates

# 2.4. Run the simulation tool executable
# =======================================

# must be exchanged for a different tool

# 2.4.1. Information
# ------------------

if [[ $verbose -gt 5 ]]; then 
    echo " gmx mdrun:  Use particle decomposition ()"
    echo "         Output log file (-g) :    " $working_dir/subdir/energy_minimization/e-min.log
    echo "         Output energy file (-e) : " $working_dir/subdir/energy_minimization/e-min.edr
    echo "         Output structure (-c) :   " $working_dir/subdir/energy_minimization/e-min.gro
    echo "         standard output:          " e-min_gmx_mdrun
    echo "" 
fi

# 2.4.2. Execution
# ----------------

   if [[ $DISTRIBUTION == "y" ]]; then
       echo gmx mdrun -nt 1 -g $working_dir/subdir/energy_minimization/e-min.log -e $working_dir/subdir/energy_minimization/e-min.edr -c $working_dir/subdir/energy_minimization/e-min.gro -s $working_dir/subdir/energy_minimization/e-min.tpr -o $working_dir/subdir/energy_minimization/e-min.trr
       `eval echo   gmx mdrun -nt 1 -g $working_dir/subdir/energy_minimization/e-min.log -e $working_dir/subdir/energy_minimization/e-min.edr -c $working_dir/subdir/energy_minimization/e-min.gro -s $working_dir/subdir/energy_minimization/e-min.tpr -o $working_dir/subdir/energy_minimization/e-min.trr` >& e-min_gmx_mdrun
   else
       gmx mdrun >& e-min_gmx_mdrun -g $working_dir/subdir/energy_minimization/e-min.log -e $working_dir/subdir/energy_minimization/e-min.edr -c $working_dir/subdir/energy_minimization/e-min.gro -s $working_dir/subdir/energy_minimization/e-min.tpr -o $working_dir/subdir/energy_minimization/e-min.trr
   fi
#fi

# 2.5. Return to the parent directory
# ===================================

cd ..

# *****************************************************************************
# 3. Pre-Pre-Equilibration
# *****************************************************************************

# The energy minimization is performed followed by a short simulation where the time step size is reduced with
# respect to the final value and the box volume is kept constant.
# The pre-pre-equilibration requires a separate ini file.
# preparation and simulation must be provided for the simulation tool

# 3.1. Main title
# ===============

if [[ $verbose -gt 3 ]]; then
    echo " Pre-pre-equilibration"
    echo " ======================"
    echo " "
fi

# 3.2. Make pre-pre-equilibration directory
# =========================================

if ! test -e pre-pre-equilibration; then
    mkdir pre-pre-equilibration
fi

cd pre-pre-equilibration


# 3.3 Update Parameters for the pre-pre-equilibration
# ===================================================

sed "s/\(^ref_t\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $temperature/" $PREPREEQUI > $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.mdp


# 3.4. Preparation routines for simulation tool
# =============================================

# must be exchanged for a different tool

#if ! test -e $working_dir/pre-pre-equi.edr || ! test -e $working_dir/pre-pre-equi.gro; then

# 3.4.1 Information
# -----------------

    if [[ $verbose -gt 7 ]]; then 
       echo " gmx grompp: Input md parameters (-f): " $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.mdp
       echo "         Input topology (-p) :     " $working_dir/current.top
       echo "         input structure (-c):     " $working_dir/subdir/energy_minimization/e-min.gro
       echo "         standard output:          " pre-pre-equi_gmx_grompp
       echo "" 
    fi


# 3.4.2. Execution
# ----------------

    gmx grompp >& pre-preequi_gmx_grompp -f $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.mdp  \
        -p $working_dir/current.top     \
        -o $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.tpr  \
        -c $working_dir/subdir/energy_minimization/e-min.gro


# 3.5. Run the simulation tool executable
# =======================================

# must be exchanged for a different tool

# 3.5.1 Information
# -----------------

if [[ $verbose -gt 5 ]]; then 
    echo " gmx mdrun:  Use particle decomposition ()"
    echo "         Output log file (-g)    : " $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.log
    echo "         Output energy file (-e) : " $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.edr
    echo "         Output structure (-c) :   " $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.gro
    echo "         standard output:          " pre-preequi_gmx_mdrun
    echo "" 
fi

# 3.5.2 Execution
# ---------------

    if [[ $DISTRIBUTION == "y" ]]; then
       echo   gmx mdrun  -nt 32 -g  $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.log -e $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.edr -c $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.gro -s $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.tpr -o $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.trr
       `eval echo   gmx mdrun  -nt 32 -g $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.log -e $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.edr -c $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.gro -s $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.tpr -o $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.trr` >& pre-preequi_gmx_mdrun
    else
       gmx mdrun >& pre-preequi_gmx_mdrun -g $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.log -e $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.edr -c $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.gro -s $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.tpr -o $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.trr
    fi
#fi

# 3.6. Return to the parent directory
# ===================================

cd ..

# *****************************************************************************
# 4. Pre-Equilibration
# *****************************************************************************
#
# (long) pre-equilibration simulations with normal time step are executed
# to shorten the equilibration run
# modification of ini file, preparation, and simulation must be provided for the simulation tool

# 4.1. Main title
# ===============

if [[ $verbose -gt 3 ]]; then
    echo " Pre-equilibriation"
    echo " =================="
    echo " "
fi
    
# 4.2. Make pre-equilibration directory
# =====================================
if ! test -e pre-equilibration; then
   mkdir pre-equilibration
fi

cd pre-equilibration

# 4.3 Update Parameters for the preequilibration
# ==============================================

# 4.3.1. Set time step, number of time steps, and restart time, if necessary
# --------------------------------------------------------------------------

# must be exchanged for a different tool

#if ! test -e $working_dir/pre-equi.edr || ! test -e $working_dir/pre-equi.gro; then
   
    p_timestep=`awk -v dt=$timestep -v ptf=$preequi_timestep_fac 'BEGIN {printf "%7.3lf\n", dt * ptf}' `
    p_nofsteps=`awk -v nsteps=$nof_steps -v psd=$preequi_steps_fac 'BEGIN {printf "%d\n", nsteps * psd}' `
    restart_time=`awk -v  dt=$timestep -v nsteps=$nof_steps -v cycle=$cycle 'BEGIN {printf "%lf\n", dt * nsteps * cycle}'`   
 
# 4.3.2. Update ini file according to commandline options
# -------------------------------------------------------

# must be exchanged for a different tool
    
    sed "s/\(^ref_t\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $temperature/" $sim_parameters |\
        sed "s/\(^gen_temp\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $temperature/" |\
        sed "s/\(^ref_p\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $pressure/" |\
        sed "s/\(^tinit\)\(.*\)\(=\)\( *\)\([0-9]*.*\)/\1\2\3 $restart_time/" |\
        sed "s/\(^nsteps\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $p_nofsteps/" |\
        sed "s/\(^nstvout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $p_nofsteps/" |\
        sed "s/\(^nstxout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $p_nofsteps/" |\
        #sed "s/\(^DispCorr\)\(.*\)\(=\)\( *\)\([A-Z,a-z]*\)/\1\2\3 no/" 
        sed "s/\(^dt \)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $p_timestep/"> $working_dir/subdir/pre-equilibration/pre-equi.mdp
    
# 4.4. Preparation routines for simulation tool
# =============================================

# must be exchanged for a different tool

# 4.4.1 Information
# -----------------

if [[ $verbose -gt 5 ]]; then 
    echo " gmx grompp: Input md parameters (-f): " $working_dir/subdir/pre-equilibration/pre-equi.mdp
    echo "         Input topology (-p) :     " $working_dir/current.top
    echo "         input structure (-c):     " $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.gro
    echo "         standard output:          " pre-equi_gmx_grompp
    echo "" 
fi

# 4.4.2. Execution
# ----------------

    gmx grompp >& preequi_gmx_grompp -f $working_dir/subdir/pre-equilibration/pre-equi.mdp \
        -p $working_dir/current.top     \
        -c $working_dir/subdir/pre-pre-equilibration/pre-pre-equi.gro  \
        -o $working_dir/subdir/pre-equilibration/pre-equi.tpr \
        
    
# 4.5. Run the simulation tool executable
# =======================================

# 4.5.1 Information
# -----------------

if [[ $verbose -gt 5 ]]; then 
    echo " gmx mdrun:  Use particle decomposition ()"
    echo "         Output log file (-g)    : " $working_dir/subdir/pre-equilibration/pre-equi.log
    echo "         Output energy file (-e) : " $working_dir/subdir/pre-equilibration/pre-equi.edr
    echo "         Output structure (-c) :   " $working_dir/subdir/pre-equilibration/pre-equi.gro
    echo "         standard output: pre-equi_gmx_mdrun "
    echo "" 
fi

# 4.5.2 Execution
# ---------------

    if [[ $DISTRIBUTION == "y" ]]; then
       `eval echo   gmx mdrun  -nt 32 -g  $working_dir/subdir/pre-equilibration/pre-equi.log -e $working_dir/subdir/pre-equilibration/pre-equi.edr -c $working_dir/subdir/pre-equilibration/pre-equi.gro -s $working_dir/subdir/pre-equilibration/pre-equi.tpr -o $working_dir/subdir/pre-equilibration/pre-equi.trr` >& preequi_gmx_mdrun
    else
       gmx mdrun >& preequi_gmx_mdrun -g  $working_dir/subdir/pre-equilibration/pre-equi.log -e $working_dir/subdir/pre-equilibration/pre-equi.edr -c $working_dir/subdir/pre-equilibration/pre-equi.gro  -s $working_dir/subdir/pre-equilibration/pre-equi.tpr -o $working_dir/subdir/pre-equilibration/pre-equi.trr
    fi

#fi

# 4.6. Return to the parent directory
# ===================================

cd ..


if [[ $equilib == "y" ]]; then

# *****************************************************************************
# 5. Equilibration cycles
# *****************************************************************************
#
# (short) equilibration simulations are executed (gmx grompp/gmx mdrun),
# evaluated and analyzed with an R script until the equilibration criterion is fulfilled
# file updating, preparation, simulation, and analysis must be modified for a different simulation tool

# 5.1. Main title
# ===============

    if [[ $verbose -gt 3 ]]; then
       echo " Equilibration cycles"
       echo " ===================="
       echo ""
    fi

# 5.2. Preparing the cycles
# =========================

# 5.2.1 Create equilibration directory and go there
# -------------------------------------------------

if ! test -e equilibration; then
    mkdir equilibration
fi

    cd equilibration

if ! test -e equi_0; then
    mkdir equi_0
fi


# 5.2.2 Redefine simulation variables
# -----------------------------------

# set energy output frequency depending on number of simulation steps
# must be exchanged for a different tool

    nstxtcout=$tmax
    nstenergy=$emax
    nstvout=$nof_steps
    nstxout=$nof_steps

# 5.2.3 Create evaluation files for simulation tool, if necessary
# ---------------------------------------------------------------

# must be exchanged for a different tool

    echo "$quantity_of_interest" > quantity.$quantity_of_interest
    if [[ $DENSITY == "y" ]]; then echo "Density" > quantity.Density; fi
    if [[ $NONBOND == "y" ]]; then echo "LJ" > quantity.LJ; fi    
    if [[ $NONBOND == "y" ]]; then echo "Coulomb" > quantity.Coulomb; fi 
    if [[ $NONBOND == "y" ]]; then echo "Coul.-recip." > quantity.Coul.-recip.; fi
    if [[ $NONBOND == "y" ]]; then echo "Disper.-corr." > quantity.Disper.-corr.; fi

# 5.2.4 data_all (Potential energy)
# ---------------------------------
    
    if ! test -e $working_dir/data_all; then
       touch $working_dir/data_all
    else
       rm $working_dir/data*
    fi

# 5.3. Start cycles
# =================

    # Stopping criterion is a well equilibrated system
    crit=FALSE
  
    # initial coordinates and energy files = output of preequilibration
    # file names must be exchanged for a different tool
    cp $working_dir/subdir/pre-equilibration/pre-equi.gro $working_dir/subdir/equilibration/equi_0/equi\_0.gro
    cp $working_dir/subdir/pre-equilibration/pre-equi.edr $working_dir/subdir/equilibration/equi_0/equi\_0.edr

    while [[  $crit != "TRUE" ]]; do

        let cycle2=cycle+1

    if ! test -e equi_$cycle2; then
       mkdir equi_$cycle2
    fi

    cd equi_$cycle2
       
        # 5.3.1 Modifying the simulation ini file
        # ---------------------------------------

        # must be exchanged for a different tool
       
        restart_time=`awk -v  dt=$timestep -v nsteps=$nof_steps -v cycle=$cycle 'BEGIN {printf "%lf\n", dt * nsteps * cycle}'`
        
        sed "s/\(^ref_t\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $temperature/" $sim_parameters |\
            sed "s/\(^gen_temp\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $temperature/" |\
            sed "s/\(^ref_p\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $pressure/" |\
            sed "s/\(^nsteps\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nof_steps/" |\
            sed "s/\(^nstenergy\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstenergy/" |\
            sed "s/\(^tinit\)\(.*\)\(=\)\( *\)\([0-9]*.*\)/\1\2\3 $restart_time/" |\
            sed "s/\(^nstvout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstvout/" |\
            sed "s/\(^nstxout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstxout/" |\
            sed "s/\(^nstxtcout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstxtcout/" |\
            sed "s/\(^dt \)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $timestep/"> $working_dir/subdir/equilibration/equi.mdp 

       if test -e $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.gro; then

          # 5.3.2 Title of cycle
          # --------------------
    
          if [[ $verbose -gt 3 ]]; then
	     echo "    Equilibration cycle " $cycle
	     echo "    ---------------------"
	     echo " "
          fi

#          if ! test -e $working_dir/equi\_$cycle.edr || ! test -e $working_dir/equi\_$cycle.gro || ! test -e $working_dir/equi\_$cycle.xtc; then


            # 5.3.3 Clean files
            # -----------------

            if test -e $working_dir/topol.tpr; then
               rm $working_dir/topol.tpr >& /dev/null
            fi

            # 5.3.4 Preparation routines for simulation tool
            # ----------------------------------------------

            # must be exchanged for a different tool
            
            # 5.3.4.1 Information
    
            if [[ $verbose -gt 5 ]]; then 
	       echo " gmx grompp: Input md parameters (-f): " $working_dir/subdir/equilibration/equi.mdp
	       echo "         Input topology (-p) :     " $working_dir/current.top
	       echo "         Input structure (-c):     " $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.gro
               echo "         Input energy (-e):        " $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.edr
	       echo "         standard output:          " equi_gmx_grompp
      	       echo "" 
            fi
    
            # 5.3.4.2 Execution

            gmx grompp >& equi_gmx_grompp -f $working_dir/subdir/equilibration/equi.mdp \
                -p $working_dir/current.top \
                -c $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.gro \
                -e $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.edr \
                -o $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.tpr 
            
            # 5.3.5 Run the simulation tool executable
            # ----------------------------------------
    
            # must be exchanged for a different tool

            # 5.3.5.1 Information
    
            if [[ $verbose -gt 5 ]]; then 
	       echo " gmx mdrun:  Use particle decomposition ()"
	       echo "         Output log file (-g)    : " $working_dir/subdir/equilibration/equi.log
	       echo "         Output energy file (-e) : " $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr
	       echo "         Output structure (-c) :   " $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.gro
               echo "         Output trajectory (-x) :  " $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.xtc
	       echo "         standard output: equi_gmx_mdrun "
	       echo "" 
            fi
    
            # 5.3.5.2 Execution 

            if [[ $DISTRIBUTION == "y" ]]; then
               `eval echo   gmx mdrun  -nt 32 -g $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.log -e $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr -c $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.gro -s $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.tpr -o $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.trr -x $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.xtc` >& equi_gmx_mdrun\_$cycle2
            else
               gmx mdrun >& equi_gmx_mdrun   -append \
                -g $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.log \
                -e $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr \
                -c $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.gro \
                -s $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.tpr \
                -o $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.trr \
                -x $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.xtc
            fi

          

            # 5.3.6 Evaluating the current cycle
            # ----------------------------------

            # all analysis routines must be exchanged for a different simulation tool
            
            # 5.3.6.1 Title
 
            if [[ $verbose -gt 5 ]]; then
	       echo "    Evaluation of equilibration cycle " $cycle2
	       echo "    -----------------------------------"
	       echo " "
            fi

            # 5.3.6.2 Get general property of interest (potential energy)

            # If you don't change anything here, the potential energy is always optimized (other properties can be implemented accordingly)
            
            if [[ $DISTRIBUTION == "y" ]]; then
               `eval echo   gmx energy -o $working_dir/$quantity_of_interest.xvg \
                -f $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr` < $working_dir/subdir/equilibration/quantity.$quantity_of_interest >& equi_gmx_energy
            else
               gmx energy >& equi_gmx_energy -o $working_dir/$quantity_of_interest.xvg \
                -f $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr < $working_dir/subdir/equilibration/quantity.$quantity_of_interest
            fi

            # 5.3.6.3 Get density

            if [[ $DENSITY == "y" ]]; then
               if [[ $DISTRIBUTION == "y" ]]; then
                  `eval echo   gmx energy -o $working_dir/Density.xvg \
                -f $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr` < $working_dir/subdir/equilibration/quantity.Density >& equi_gmx_energy_density
               else
                  gmx energy >& equi_gmx_energy_density -o $working_dir/Density.xvg \
                -f $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr < $working_dir/subdir/equilibration/quantity.Density
               fi
            fi

            # 5.3.6.4 Get nonbonded energy

            if [[ $NONBOND == "y" ]]; then
               $MOSCITOBINPATH/trajenergy -sys $working_dir/current.sys -par $MOSCITOPARFILE -in $diff_structure \
               -xtc $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.xtc > $working_dir/traj_energy_equi\_$cycle2.dat     

               grep EPOT_INTER_SUM $working_dir/traj_energy_equi\_$cycle2.dat | awk '{print $2}' > $working_dir/Intermolecular\_$cycle2.list
            fi
            
            # 5.3.6.5 Convert files into GROW standard lists

            # must be modified, if current output format differs from GROW standard lists            

            egrep -v "\@|#" $working_dir/$quantity_of_interest.xvg > $working_dir/$quantity_of_interest\_$cycle2.list
            
            if [[ $DENSITY == "y" ]]; then
               egrep -v "\@|#" $working_dir/Density.xvg > $working_dir/Density\_$cycle2.list;
            fi
            
            # 5.3.7 Analyse the current data
            # ------------------------------

            # 5.3.7.1 Title

            if [[ $verbose -gt 5 ]]; then
	       echo "    Analysis of equilibration cycle " $cycle2
	       echo "    --------------------------------"
	       echo " "
            fi     
            
            # 5.3.7.2 Append the evaluated data to a total equilibration evaluation file
            
            cat  $working_dir/$quantity_of_interest\_$cycle2.list >>  $working_dir/data_all   
            rm  $working_dir/$quantity_of_interest\_$cycle2.list >& /dev/null
 
            if [[ $DENSITY == "y" ]]; then
               cat  $working_dir/Density\_$cycle2.list >>  $working_dir/data_all_density;
               rm  $working_dir/Density\_$cycle2.list >& /dev/null;
            fi
            if [[ $NONBOND == "y" ]]; then
               cat  $working_dir/Intermolecular\_$cycle2.list >>  $working_dir/data_all_nonbond;   
               rm  $working_dir/Intermolecular\_$cycle2.list >& /dev/null;
            fi
            
            # 5.3.7.3 Extract from all data the last $nof_samples data points            

            all_data=`wc -l $working_dir/data_all | awk '{print $1}'`
            
            tail -$nof_samples $working_dir/data_all > $working_dir/data

            if [[ $DENSITY == "y" ]]; then
               tail -$nof_samples $working_dir/data_all_density > $working_dir/data_density;
            fi
            if [[ $NONBOND == "y" ]]; then
               tail -$nof_samples $working_dir/data_all_nonbond > $working_dir/data_nonbond;
            fi     

            # 5.3.7.4 Do the statistics via R

            if [[ $verbose -gt 3 ]]; then  
	       echo " Statistics via R "
	       echo "" 
            fi

            crit=`R  --slave < $r_script --args $working_dir/data_all $cycle2 Potential $SPOT $TW | awk '{print $NF}'`
                 
            if [[ $DENSITY == "y" ]]; then
               crit_density=`R  --slave < $r_script --args $working_dir/data_all_density $cycle2 Density $SRHO $TW | awk '{print $NF}'`;
            else
               crit_density=TRUE;
            fi

            if [[ $NONBOND == "y" ]]; then
               crit_nonbond=`R  --slave < $r_script --args $working_dir/data_all_nonbond $cycle2 Nonbond $SNB $TW | awk '{print $NF}'`;
            else
               crit_nonbond=TRUE;
            fi

            if [[ $crit == "TRUE" && $crit_density == "TRUE" && $crit_nonbond == "TRUE" ]]; then
               crit=TRUE;
            else
               cycle=`expr $cycle + 1`
#               cp $working_dir/subdir/equilibration/equi\_$cycle.gro $working_dir/subdir/equilibration/equi\_$cycle2.gro
               crit=FALSE;
            fi

            echo "             Equilibration criterion fulfilled? $crit"
            
            rm \#* >& /dev/null 
            
            if test -e $working_dir/step*.pdb; then
                
                echo "DANGER step file"
                exit='exit'
            fi
                 
        else
            echo "file $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.gro missing"
            echo "file $working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.gro missing" > $working_dir/fatal_error
            exit 1
        fi
        
        
        if test '$exit' == 'exit'; then
            exit
        fi

        cd ..
        
    done
    
    if test '$exit' == 'exit'; then
        exit
    fi

    # 5.3.7.5 Finishing equilibration: writing final coordinates file for production
    
    cp $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.gro $working_dir/system_equilibrated.gro
    cp $working_dir/subdir/equilibration/equi_$cycle2/equi\_$cycle2.edr $working_dir/system_equilibrated.edr
    
    if [[ $verbose -gt 3 ]]; then
       echo "equilibration finished: $crit"
    fi

else

    cp $working_dir/subdir/pre-equilibration/pre-equi.gro $working_dir/system_equilibrated.gro
    cp $working_dir/subdir/pre-equilibration/pre-equi.edr $working_dir/subdir/equilibration/equi\_$cycle2.edr


fi # end of equilibration

# 5.4. Return to the parent directory
# ===================================

cd ..

fi

# *****************************************************************************
# 6. Production run
# *****************************************************************************

# 6.1. Main title
# ===============

if [[ $verbose -gt 3 ]]; then
    echo " Production run"
    echo " =============="
    echo ""
fi

# 6.2. Make production directory
# ==============================

if ! test -e production; then
   mkdir production
fi

cd production


#  if ! test -e $working_dir/production/production.edr || ! test -e $working_dir/production/production.gro || ! test -e $working_dir/production/production.xtc; then

# 6.3. Modifying simulation ini file
# =================================

# 6.3.1. Set time step, number of time steps, and restart time, if necessary
# --------------------------------------------------------------------------

# must be exchanged for a different tool

    
    nof_steps=`awk -v nsteps=$nof_steps -v dsf=$prod_steps_fac 'BEGIN {printf "%d\n", nsteps * dsf}'`
    nstxtcout=$tmax
    nstenergy=$emax
    nstvout=$nof_steps
    nstxout=$nof_steps
 
# 6.3.2. Update ini file according to commandline options
# -------------------------------------------------------

# must be exchanged for a different tool

    if [[ $RERUN == "y" ]]; then
       nof_steps=`awk -v nsteps=$nof_steps -v freq=$tmax 'BEGIN {printf "%d\n", nsteps/freq}'`
       nstenergy=1
       nstvout=0
       nstxout=0
       nstxtcout=0
    fi

    sed "s/\(^ref_t\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $temperature/" $sim_parameters |\
        sed "s/\(^ref_p\)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $pressure/"  |\
        sed "s/\(^nsteps\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nof_steps/"  |\
        sed "s/\(^nstenergy\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstenergy/"   |\
        sed "s/\(^nstvout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstvout/" |\
        sed "s/\(^nstxout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstxout/" |\
        sed "s/\(^nstxtcout\)\(.*\)\(=\)\( *\)\([0-9]*\)/\1\2\3 $nstxtcout/" |\
        sed "s/\(^dt \)\(.*\)\(=\)\( *\)\([0-9\.]*\)/\1\2\3 $timestep/"> $working_dir/subdir/production/production.mdp
 
# 6.4. Preparation routines for simulation tool
# =============================================

# must be exchanged for a different tool

# 6.4.1 Information
# -----------------

if [[ $rerun_exec == "TRUE" ]]; then
       loop=`awk -v l=$loop 'BEGIN {printf "%d\n", l-1}'`
       gro_file=$TMPPATH/$TEMPERATURE/original.out.$loop/subdir/production/production.gro
       edr_file=$TMPPATH/$TEMPERATURE/original.out.$loop/subdir/production/production.edr
else
       gro_file=$working_dir/system_equilibrated.gro
       edr_file=$working_dir/subdir/equilibration/equi_$cycle/equi\_$cycle.edr
fi

if [[ $verbose -gt 5 ]]; then 
    echo " gmx grompp: Input md parameters (-f): " $working_dir/subdir/production/production.mdp
    echo "         Input topology (-p) :     " $working_dir/current.top
    echo "         Input structure (-c):     " $gro_file
    echo "         Input energy (-e):        " $edr_file
    echo "         standard output:          " prod_gmx_grompp
    echo "" 
fi

# 6.4.2. Execution
# ----------------    

    

    gmx grompp >& prod_gmx_grompp -f $working_dir/subdir/production/production.mdp -p $working_dir/current.top  \
        -c $gro_file \
        -e $edr_file \
        -o $working_dir/subdir/production/production.tpr  

# 6.5. Run the simulation tool executable
# =============================================

# must be exchanged for a different tool

# 6.5.1. Information
# ------------------

if [[ $verbose -gt 5 ]]; then
    echo " gmx mdrun:  Use particle decomposition ()"
    echo "         Output log file (-g) :    " $working_dir/subdir/production/production.log
    echo "         Output energy file (-e) : " $working_dir/subdir/production/production.edr
    echo "         Output structure (-c) :   " $working_dir/subdir/production/production.gro
    echo "         standard output:          " prod_gmx_mdrun
    echo ""
fi

# 6.5.2. Execution
# ----------------

            if [[ $rerun_exec == "TRUE" ]]; then
               traj_file=$TMPPATH/$TEMPERATURE/original.out.$loop/subdir/production/production.xtc
               rerun_option="-rerun"
            fi

            if [[ $DISTRIBUTION == "y" ]]; then
               `eval echo   gmx mdrun  -nt 32   -append $rerun_option $traj_file \
                -g $working_dir/subdir/production/production.log \
                -e $working_dir/subdir/production/production.edr \
                -c $working_dir/subdir/production/production.gro \
                -s $working_dir/subdir/production/production.tpr \
                -x $working_dir/subdir/production/production.xtc` >& prod_gmx_mdrun
            else
               gmx mdrun >& equi_gmx_mdrun   -append $rerun_option $traj_file \
                -g $working_dir/subdir/production/production.log \
                -e $working_dir/subdir/production/production.edr \
                -c $working_dir/subdir/production/production.gro \
                -s $working_dir/subdir/production/production.tpr \
                -x $working_dir/subdir/production/production.xtc
            fi

            if [[ $rerun_exec == "FALSE" ]]; then
               traj_file=$working_dir/subdir/production/production.xtc
            fi

# 6.6. Evaluate production run
# =============================================

    
# 6.6.1 Calculate density, temperature, volume, pressure, and potential energy
# ----------------------------------------------------------------------------

# must be exchanged for a different tool


for quantity_of_interest in LJ Coulomb Disper.-corr. Coul.-recip. Density Pressure Temperature Volume Potential; do
        
        echo "$quantity_of_interest" > quantity.$quantity_of_interest 
        
        if [[ $DISTRIBUTION == "y" ]]; then
           `eval echo   gmx energy -o $working_dir/$quantity_of_interest.xvg \
            -f $working_dir/subdir/production/production.edr` < quantity.$quantity_of_interest >& prod_gmx_energy
        else
           gmx energy >& prod_gmx_energy -o $working_dir/$quantity_of_interest.xvg \
            -f $working_dir/subdir/production/production.edr < quantity.$quantity_of_interest
        fi
        
        egrep -v "\@|#" $working_dir/$quantity_of_interest.xvg > $working_dir/$quantity_of_interest.list
        
        
    done

cp $working_dir/Density.list $working_dir/Bulk_Density.list
cat $working_dir/Volume.list | awk '{print $2/10**27}' > $working_dir/Volume2.list
mv $working_dir/Volume2.list $working_dir/Volume.list



# 6.6.2 Calculate diffusion coefficient
# ----------------------------------------------------------------------------


if [[ $pressure == 1.00001 ]]; then

        time=`echo "$tmax*0.002" | bc`
        msd=`echo "$diff_tmax/4" | bc`

        if [[ $DISTRIBUTION == "y" ]]; then
            echo "  /home/mhuelsma/software/moscito-4.150/bin/msdmol"
           `eval echo   /home/mhuelsma/software/moscito-4.150/bin/msdmol  \
            -sys $working_dir/current.sys -gap 2 \
            -in $diff_structure \
            -time $time \
            -tmax $diff_tmax \
            -xtc $traj_file` >& $working_dir/subdir/production/prod_msd
        else
           ~/software/moscito-4.150/bin/msdmol >& prod_moscito  \
            -sys $working_dir/current.sys \
            -in $diff_structure \
            -time $time \
            -tmax $diff_tmax \
            -xtc $traj_file
        fi


           tail -$msd msd.dat | awk '{print $1,$2}' | $FITPOLYBINPATH/fitpoly -o 1 | tail -4 | head -1 | awk '{print $3/6*1000}' > $working_dir/Density.list


     
fi

# must be exchanged for a different tool

    if [[ $DIFFCOEF == "y" ]]; then

        time=`echo "$tmax*0.002" | bc`
        msd=`echo "$diff_tmax/4" | bc`

        if [[ $DISTRIBUTION == "y" ]]; then
           `eval echo   $MOSCITOBINPATH/msdmol  \
            -sys $working_dir/current.sys \
            -in $diff_structure \
            -time $time \
            -tmax $diff_tmax \
            -xtc $traj_file` >& $working_dir/subdir/production/prod_msd
        else
           $MOSCITOBINPATH/msdmol >& prod_moscito  \
            -sys $working_dir/current.sys \
            -in $diff_structure \
            -time $time \
            -tmax $diff_tmax \
            -xtc $traj_file
        fi

        if [[ $DIFFCOEFCAT == "y" ]]; then

           tail -$msd msd.dat | awk '{print $1,$2}' | $FITPOLYBINPATH/fitpoly -o 1 | tail -4 | head -1 | awk '{print $3/6*1000}' > $working_dir/D_cat.erg

        fi

        if [[ $DIFFCOEFAN == "y" ]]; then

           tail -$msd msd.dat | awk '{print $1,$3}' | $FITPOLYBINPATH/fitpoly -o 1 | tail -4 | head -1 | awk '{print $3/6*1000}' > $working_dir/D_an.erg 

        fi

     fi


# 6.6.3 Calculate reorientation time
# ----------------------------------

# must be exchanged for a different tool


     if [[ $REOR == "y" ]]; then

	reor_tmax=`echo "$diff_tmax*2" | bc`

        time=`echo "$tmax*0.002" | bc`
        $MOSCITOBINPATH/vectorcor >& $working_dir/subdir/production/prod_vectorcor -sys $working_dir/current.sys -in $diff_structure -xtc $traj_file -time $time -tmax $reor_tmax

        laeng=`echo "$reor_tmax/4" | bc`  

        laeng2=`echo "$laeng*3/4" | bc`  

        tau2=`cat $working_dir/subdir/production/vec_p1p2.dat | awk '{print $1,$3}' | $MOSCITOBINPATH/integrate | awk '{print $2}' | tail -$laeng | head -$laeng2  | $MOSCITOBINPATH/average`

        echo $tau2 | awk '{print $1}' | $MOSCITOBINPATH/average | awk '{print $1}' > $working_dir/tau2.erg
     fi



# 6.6.4 Calculate viscosity
# ----------------------------------

# must be exchanged for a different tool

     if [[ $VISC == "y" ]]; then

 if test -e visc.dat; then
    rm visc.dat
 fi
 time=`echo "$emax*0.002" | bc`
 nvalue=`echo "$nof_steps/$emax" | bc`
 nzero=`echo "$nvalue-1000" | bc`
 cat=`echo "1000/1.2" | bc`
 cat2=`echo "$cat*3/4" | bc`


cat > inp <<EOF
'p_tensor.list'    input
'p0pt'             output
$nvalue            nvalue
$nzero             nzero
$time              deltat
EOF

echo 27 28 29 30 31 32 33 34 35 0 | gmx energy >& prod_visc -f $working_dir/subdir/production/production.edr -o $working_dir/subdir/production/energy.xvg
tail -$nvalue $working_dir/subdir/production/energy.xvg | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' > $working_dir/subdir/production/p_tensor.list
rm $working_dir/subdir/production/energy.xvg


$MOSCITOBINPATH/viscosity_xy < $working_dir/subdir/production/inp

 cat $working_dir/subdir/production/p0pt | awk '{if (NR == 1) {y = -$2*0.01; x = 0.02} y += $2*0.02;printf "%lf %.7lf\n", $1+x, y}' > $working_dir/subdir/production/p0pt_xy.int.dat
# tail -$cat $working_dir/subdir/production/p0pt_xy.int.dat | head -$cat2 >> $working_dir/subdir/production/visc.tmp
 tail -$cat $working_dir/subdir/production/p0pt_xy.int.dat | head -$cat2 | awk '{print $2}' | $MOSCITOBINPATH/average | awk '{print $1}' >> $working_dir/subdir/production/visc.dat

$MOSCITOBINPATH/viscosity_xz < $working_dir/subdir/production/inp

 cat $working_dir/subdir/production/p0pt | awk '{if (NR == 1) {y = -$2*0.01; x = 0.02} y += $2*0.02;printf "%lf %.7lf\n", $1+x, y}' > $working_dir/subdir/production/p0pt_xz.int.dat

 tail -$cat $working_dir/subdir/production/p0pt_xz.int.dat | head -$cat2 | awk '{print $2}' | $MOSCITOBINPATH/average | awk '{print $1}' >> $working_dir/subdir/production/visc.dat

$MOSCITOBINPATH/viscosity_yz < $working_dir/subdir/production/inp

 cat $working_dir/subdir/production/p0pt | awk '{if (NR == 1) {y = -$2*0.01; x = 0.02} y += $2*0.02;printf "%lf %.7lf\n", $1+x, y}' > $working_dir/subdir/production/p0pt_yz.int.dat

 tail -$cat $working_dir/subdir/production/p0pt_yz.int.dat | head -$cat2 | awk '{print $2}' | $MOSCITOBINPATH/average | awk '{print $1}' >> $working_dir/subdir/production/visc.dat

 R --slave < $VISCSCRIPT --args $working_dir/subdir/production/visc.dat $working_dir/Volume.list $working_dir/Temperature.list $working_dir


#rm $working_dir/subdir/production/p_tensor.list $working_dir/subdir/production/p0pt $working_dir/subdir/production/inp $working_dir/subdir/production/visc.dat
   fi

# 6.6.5 Calculate nonbonded energy and enthalpy of vaporization
# -------------------------------------------------------------

# must be exchanged for a different tool

if [[ $NONBOND == "y" ]]; then

   $MOSCITOBINPATH/trajenergy -sys $working_dir/current.sys -par $MOSCITOPARFILE -in $diff_structure \
   -xtc $traj_file > $working_dir/traj_energy.dat     

   grep EPOT_INTER_SUM $working_dir/traj_energy.dat | awk '{print $2}' > $working_dir/Intermolecular.list

   grep EPOT_INTRA_SUM $working_dir/traj_energy.dat | awk '{print $2}' > $working_dir/Intramolecular.list 
   sleep 60
   R --slave < $VAPORSCRIPT --args $working_dir/Intermolecular.list $NMOL $working_dir/Temperature.list $working_dir $CORRECTION

fi

# 6.6.6 Calculate molecular weight
# --------------------------------

# must be exchanged for a different tool


mass_of_system=`awk -v m=$MASS -v N=$NMOL 'BEGIN {printf "%7.50lf\n", m*N/(6.022*10**26)}' `
echo $mass_of_system > $working_dir/mass_of_system.list
        
       

# 6.6.7 Calculate other properties required for optimization workflow
# -------------------------------------------------------------------

# must be exchanged for a different tool

# *****************************************************************************
# 7. Clean up and Termination
# *****************************************************************************

rm $OUTPATH/subdir/*/*.tpr
rm $OUTPATH/subdir/*/*.trr
#rm $OUTPATH/subdir/*/*.xtc

touch $TMPPATH/*
touch $TMPPATH/*/*
touch $TMPPATH/*/*/*
touch $TMPPATH/*/*/*/*
touch $TMPPATH/*/*/*/*/*
touch $TMPPATH/*/*/*/*/*/*

echo "mkjobequigromacs.sh erfolgreich beendet"

echo "Simulation erfolgreich beendet" > $OUTPATH/terminated.txt

cd $act_dir
