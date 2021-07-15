#     abstract class Molecular_Simulation_Optimization_Problem, 2014/05/06 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		molecular_simulation_optimization_problem.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		06-05-2014
#
#     Description:	abstract class Molecular_Simulation_Optimization_Problem
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
#			
#     Output:	        
#
#     Imported:		loss.py
#			optimization_problem.py
#
#     Calling:	
#			
#
#     Modifications:
#
#	$Log: md_opt$
#
#     end_doc

# ***********************************************************************
# required python modules
# ***********************************************************************

import generic_optimization_problem.loss
from generic_optimization_problem.optimization_problem import Optimization_Problem

# ***********************************************************************
# class Molecular_Simulation_Optimization_Problem
# ***********************************************************************

class Molecular_Simulation_Optimization_Problem(Optimization_Problem):
       
    def __init__(self, dimension, physical_properties_loss, force_field_constraints):
          """ constructor """

          Optimization_Problem.__init__(self, dimension, physical_properties_loss, force_field_constraints)

    def __del__(self):
          """
          destructor
          """	
	  del self	
