#     abstract class Optimization_Problem, 2014/04/30 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		optimization_problem.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		30-04-2014
#
#     Description:	abstract class Optimization_Problem
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
#			
#     Output:	        
#
#     Imported:		box_constraints.py
#			objective_function.py
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

import box_constraints
import objective_function

# ***********************************************************************
# class Optimization_Problem
# ***********************************************************************

class Optimization_Problem:
       
    def __init__(self, dimension, objective_function, box_constraints):
          """ constructor """
          
          self._dimension = dimension
	  self._objective_function = objective_function
          self._box_constraints = box_constraints

    def get_dimension(self):
        """ return the dimension of the optimization problem """
        return self._dimension

    def get_objective_function(self):
        """ return the objective function of the optimization problem """
        #print "objective function called in",__name__
        #print "what is in self._objective_function: ",self._objective_function
        return self._objective_function

    def get_box_constraints(self):
        """ return the box constraints of the optimization problem """
        return self._box_constraints

    def is_optimal(self, x, tolerance, Fxmin = None):
	""" stopping criterion fulfilled at x? """

        if Fxmin != None:
           return Fxmin <= tolerance

        else:
	   try:
              return abs(self._objective_function.get_function_value(x)) <= tolerance
           except:
              return False

    def __del__(self):
          """
          destructor
          """	
	  del self	
