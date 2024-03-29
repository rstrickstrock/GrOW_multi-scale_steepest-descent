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

- copy the "Code" directory to the desired location
- navigate to Code/GrOW/
- exectute the main script the following way:
    python main.py
    e.g.
    $ python main.py ../opt_1/octane_hybrid.cfg

    
    
# Adapting the file to fit the IT infrastructure
    
- make sure all the relative filepaths in the -file are correct
- adapt files in
    - Code/GrOW/parallel_jobs/*
    - Code/GrOW/simulation/*
      
      to fit to your simulation tools and environments used (e.g. cluster queueuing software)
