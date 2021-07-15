#     abstract class Box_Constraints, 2014/04/30 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		box_constraints.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		30-04-2014
#
#     Description:	abstract class Box_Constraints
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
#			
#     Output:	        
#
#     Imported:		opt_class.py
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

import opt_class

opt = opt_class.Optimization()

# ***********************************************************************
# class Box_Constraints
# ***********************************************************************

class Box_Constraints:
       
    def __init__(self, dimension, initial_guess):
          """ constructor """

          if len(initial_guess) != dimension:
             print "Error in box_constraints: The length of the initial_guess (%d) does not correspond to the dimension (%d)!" %(len(initial_guess),dimension)
             opt.errorexit()


          self._dimension = dimension
          self._initial_guess = initial_guess
          
          self._min_vector = []
          self._max_vector = []

    def get_boundary(self):
        """ return the boundary vectors """
        return [self._min_vector, self._max_vector]

    def is_feasible(self, x):
        """ decides if a vector x is feasible in the sense of the box constraints """

        for i in range(len(x)):
            if x[i] < self._min_vector[i] or x[i] > self._max_vector[i]:
               return False

        return True

    def __del__(self):
          """
          destructor
          """	
	  del self	
