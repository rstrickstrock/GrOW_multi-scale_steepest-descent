#     abstract class Armijo_Step_Length_Control, 2014/05/07 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		armijo_step_length_control.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		07-05-2014
#
#     Description:	class Armijo_Step_Length_Control
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
# 			step_length_control.py
#
#     Calling:	
# 			
#
#     Modifications:
#
# 	$Log: md_opt$
#
#     end_doc

# ***********************************************************************
# required python modules
# ***********************************************************************

import generic_optimization_problem.optimization_problem
import opt_class
import os

from step_length_control import Step_Length_Control

opt = opt_class.Optimization()

# ***********************************************************************
# class Armijo_Step_Length_Control
# ***********************************************************************

class Armijo_Step_Length_Control(Step_Length_Control):
       
    def __init__(self, optimization_problem, config):
        """ constructor """
           
        Step_Length_Control.__init__(self, optimization_problem, config)

        # sigma
	try:
	   	sigma = float(config.get("OPT", "sigma"))
        except ValueError:
                print "ERROR in armijo_step_length_control: The armijo parameter sigma must be a numerical value. %s is not numerical." % (config.get("OPT", "sigma"))
		opt.errorexit()
        except:
		print "ERROR in armijo_step_length_control: You must indicate the armijo parameter sigma."
		opt.errorexit()

   	if sigma <= 0 or sigma >= 1:
      	   print "ERROR in armijo_step_length_control: The armijo parameter sigma must be in the interval (0,1). %s is not." % (sigma)
      	   opt.errorexit()

        # beta
        if os.getenv("METHOD") == "trust_region" and os.getenv("TR_partial") == "exact":
               try:
     	        	beta = float(config.get("OPT", "beta"))
	       except ValueError:
                        print "ERROR in armijo_step_length_control: The armijo parameter beta must be a numerical value. %s is not numerical." % (config.get("OPT", "beta"))
 		        opt.errorexit()
     	       except:
		        print "ERROR in armijo_step_length_control: You must indicate the armijo parameter beta."
		        opt.errorexit()

    	       if beta <= 0 or beta >= 1:
      	          print "ERROR in armijo_step_length_control: The armijo parameter beta must be in the interval (0,1). %s is not." % (beta)
                  opt.errorexit()

               if sigma >= 0.5:
                  print "ERROR in armijo_step_length_control: The armijo parameter sigma must be in the interval (0,0.5). %s is not." % (sigma)
      	    	  opt.errorexit()
                    
        # number of parallel armijo evaluations
	try:
		parallel = int(config.get("OPT", "number_of_parallel_armijo_evaluations"))
        except ValueError:
	        print "ERROR in armijo_step_length_control: The number of parallel armijo evaluations must be an integer. %s is not." % (config.get("OPT", "number_of_parallel_armijo_evaluations"))
 		opt.errorexit()
     	except:
		print "ERROR in armijo_step_length_control: You must indicate the number of parallel armijo evaluations."
		opt.errorexit()

        if parallel <= 0:
           print "ERROR in armijo_step_length_control: Your indiacted number of parallel armijo evaluations must be greater than 0."
           opt.errorexit()

        if parallel > self._number_of_iterations:
           print "ERROR in armijo_step_length_control: Your indiacted number of parallel armijo evaluations is greater than the total number of armijo iterations."
           opt.errorexit()

        self._sigma = sigma
        self._parallel = parallel

    def get_step_length(self, min_vector, max_vector, x, f_x, gradient, descent):
        """ defines an Armijo step length """

        print "Computing Armijo step length"

        # read necessary environment variables
        method = os.getenv("METHOD")
 
        # initializations
        self.feasible_step_length(min_vector, max_vector, x, descent)
        beta = self._step_length
        gradient_descent_prod = opt.scalar_product(gradient, descent)
	
        # for higher step lengths
        if beta >= 1.0 and gradient_descent_prod < 0:
           c = 0.01 - pow(opt.norm(descent), 2) * beta / opt.scalar_product(gradient, descent)
           s = -c * gradient_descent_prod / pow(opt.norm(descent), 2)
           beta = 0.5
        else:
           s = 1.0
           if beta >= 1.0:
              beta = 0.9

        # Armijo step length control
        self._number_of_iteration = 2
        self._number_of_iterations += self._number_of_iteration
        break_ind = False
        while True:
              if self._number_of_iteration == self._number_of_iterations:
                 print "Number of Armijo steps exceeded."
                 opt.errorexit()
              if os.getenv("DERIV_NORM") == None:
                 factor = 1./opt.norm(descent)
              else:
                 factor = 1.0
              parameter_set = []
              for k in range(self._parallel):
                  print "Armijo Iteration %d" % (self._number_of_iteration + k)                 
                  # new candidate iteration        
                  x_new = []
	          for i in range(len(descent)):                 
                      x_new.append(x[i] + s * pow(beta, self._number_of_iteration + k) * descent[i] * factor)
                  parameter_set.append(x_new)
              f_x_armijo_list = self._optimization_problem.get_objective_function().get_function_values(parameter_set, self._number_of_iteration)
              # armijo condition evaluation
              for f_x_armijo in f_x_armijo_list:
                  if f_x_armijo > f_x + self._sigma * s * pow(beta, self._number_of_iteration) * gradient_descent_prod * factor:
                     armijo_evaluation = "no"
                  else:
                     armijo_evaluation = "ok"
                  if armijo_evaluation == "ok" or s * pow(beta, self._number_of_iteration) < 0.0:    
                     self._step_length = s * pow(beta, self._number_of_iteration)
                     print "Computed Armijo step length t = %f" % (self._step_length)
                     if armijo_evaluation != "ok":
                        print "!!! Step length may be not efficient!!!"
                     break_ind = True
                     break
                  self._number_of_iteration += 1
              if break_ind:
                 break
        
        return [self._step_length, f_x_armijo]

    def __del__(self):
          """
          destructor
          """	
	  del self	
