#     abstract class Objective_Function, 2014/04/30 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		objective_function.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		30-04-2014
#
#     Description:	abstract class Objective_Function
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
#			
#     Output:	        
#
#     Imported:		
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
# class Objective_Function
# ***********************************************************************

class Objective_Function:
       
    def __init__(self, dimension):
          """ constructor """
          
          self._dimension = dimension

    def get_dimension(self):
        """ return the dimension of the objective function """
        return self._dimension

    # virtual methods
    
    def get_function_values(self, parameter_set):
        """ returns a list of function values for a given parameter set"""
        print "got function values in",__name__
        pass
    
    def get_function_value(self, x):
        """ return the objective function value of a vector x """
        pass

    def get_gradient(self, x):
        """ return the gradient of the objective function at x """
        pass

    def get_hessian(self, x):
        """ return the Hessian of the objective function at x """
        pass

    def __del__(self):
          """
          destructor
          """	
	  del self	
