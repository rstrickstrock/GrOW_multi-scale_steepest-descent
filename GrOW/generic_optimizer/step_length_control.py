#     abstract class Step_Length_Control, 2014/04/30 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		step_length_control.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		30-04-2014
#
#     Description:	class Step_Length_Control
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
#			
#     Output:	        
#
#     Imported:		optimization_problem.py
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

import generic_optimization_problem.optimization_problem
import opt_class
import os

opt = opt_class.Optimization()

# ***********************************************************************
# class Step_Length_Control
# ***********************************************************************

class Step_Length_Control:
       
    def __init__(self, optimization_problem, config):
          """ constructor """
           
	  try:
  	     sl_method = config.get("OPT","sl_method")
          except:
	     print "ERROR in step_length_control: You must indicate the step length method."
  	     opt.errorexit()

          if sl_method not in ["armijo"]:
             print "ERROR in step_length_control: At the moment, 'armijo' is the only valid step length control. '%s' is not valid." %(sl_method)
             opt.errorexit()

          try:
		number_of_iterations = int(config.get("OPT","max_armijo"))
          except ValueError:
                print "ERROR in step_length_control: The variable 'max_armijo' must be an integer. %s is not." %(config.get("OPT","max_armijo"))
                opt.errorexit()
          except:
		print "ERROR in step_length_control: You must indicate the maximal number of step length control iterations."
                opt.errorexit()

          self._optimization_problem = optimization_problem
          self._sl_method = sl_method
          self._number_of_iterations = number_of_iterations
	  self._step_length = None # initialization of step length
          self._number_of_iteration = 0

    def get_step_length(self):
        """ returns the actual step length """
        return self._step_length

    def set_step_length(self):
        """ virtual function for finding a step length """
        pass

    def get_number_of_iteration(self):
        """ returns the actual iteration number """
        return self._number_of_iteration

    def feasible_step_length(self, min_vector, max_vector, x, descent):
        """ computes a step length that is feasible in the sense of lowerbound <= ||x + sl * descent|| <= upperbound """

        if os.getenv("DERIV_NORM") == None:
           descent = opt.mult_vector(1/opt.norm(descent),descent)

        if descent[0] > 0:
           t = (max_vector[0]-x[0])/descent[0]
        elif descent[0] < 0:
           t = (min_vector[0]-x[0])/descent[0]
        else:
           t = 1.0

        for j in range(len(x)):
            x_new_j = x[j] + t*descent[j]
            if x_new_j < min_vector[j] or x_new_j > max_vector[j]: 
               if descent[j] > 0.0:
                  t = (max_vector[j] - x[j])/descent[j]
               if descent[j] < 0.0:
                  t = (min_vector[j] - x[j])/descent[j]

        if t == 0.0:
           print "ERROR in step_length_control: Computed step length is zero. The border of the admissible domain is reached and the direction leads out of it."
           print "No minimum found in admissible domain."
           opt.errorexit()

	self._step_length = t

    def __del__(self):
          """
          destructor
          """	
	  del self	
