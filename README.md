# GrOW_multi-scale_steepest-descent
Python and input files for reproduction of the results of COMPHY-D-21-00093R2

This repository contains program and input files for the reproduction of the results presented in COMPHY-D-21-00093R2.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    

To reproduce the results and use the given input files installations of Gromacs <https://www.gromacs.org/Gromacs_papers> and Amber <https://ambermd.org/> are necessary.

Furthermore, adaptions to the program files may be necessary depending on the used IT infrastructure. 


# Installation and execution guide

- copy the code directory ("FFLOW/") to the desired location
- navigate to FFLOW/
- exectute the main script the following way:

    python main.py <config>
    
    e.g.
    
    $ python main.py ../opt_1/octane_hybrid.cfg
  
- The optimization in "input_data/" are started via scripts (e.g. "input_data/wf_exp_init/batch_run_PP-15_MM-85.sh") that also utlize the provided IT structure at the time the optimizations were performed. 
    
    
# Adapting the file to fit the IT infrastructure
    
- adapt the paths for 'cwd' in all PP-*/octane_hybrid_new.cfg files
- adapt the paths for 'cwd' in all run_PP*.sh files
- adapt the configurations in all PP-*/batch_run_PP*.sh files to match your IT infrastructure
- maybe adapt files in 
 
        - Code/GrOW/parallel_jobs/*
		
        - Code/GrOW/simulation/*
    
  	    to fit to your simulation tools and environments used (e.g. cluster queueuing software)
