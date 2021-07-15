#     abstract class Optimization_Algorithm, 2014/04/30 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		optimization_algorithm.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		30-04-2014
#
#     Description:	abstract class Optimization_Algorithm
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
#                       optimization_problem.py
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

import os

import generic_optimization_problem.optimization_problem
import step_length_control

import opt_class
opt = opt_class.Optimization()

# ***********************************************************************
# class Optimization_Algorithm
# ***********************************************************************

class Optimization_Algorithm:
       
    def __init__(self, optimization_problem, step_length_control, number_of_iterations, config):
        
        
        """ constructor """
        self._optimization_problem = optimization_problem
        self._step_length_control = step_length_control
        self.__number_of_iterations = number_of_iterations
          
        self._x = None  # iteration
        self._number_of_iteration = 0  # actual iteration

	    # starting point
        try: 
            par = config.get("OPT", "parameters")
        except:
            print "ERROR in optimization_algortihm.py: You must indicate a parameter file."
            opt.errorexit()

	    try:
	        os.stat(par)
	    except:
	        print "ERROR in optimization_algortihm.py: Parameter file %s does not exist." % (par)
            opt.errorexit()

        self._x0 = opt.read_last_parameter(par) 

        if len(self._x0) != optimization_problem.get_dimension():
            print "ERROR in optimization_algortihm.py: The length of the parameter vector (%d) is not equal to the dimension of the optimization problem (%d)." % (len(self._x0), optimization_problem.get_dimension())
            opt.errorexit()

        self.set_iteration(self._x0)

	    # tolerance level
        try:
            limit = float(config.get("OPT", "limit"))
        except ValueError:
            print "ERROR in optimization_algortihm.py: The tolerance level '%s' is not a numerical value." % (config.get("OPT", "limit"))
       	    opt.errorexit()
        except:
            print "ERROR in optimization_algortihm.py: You must indicate the tolerance level."
            opt.errorexit()

        opt.setenv("LIMIT", repr(limit)) 
        self.__tolerance = limit 

        # function evaluation only or complete optimization?
        try:
            self.__workflow = config.get("OPT", "workflow")
        except:
            print "ERROR in optimization_algorithm.py: You must indicated the worklow, i.e. whether you want to perform an optimization or one function evaluation only."
            opt.errorexit()

        if self.__workflow not in ["optimization", "feval"]:
            print "ERROR in optimization_algorithm.py: For the workflow you must decide between 'optimization' and 'feval'. '%s' is not valid." % (self.__workflow)
            opt.errorexit()

        try:
            self._objective_function = config.get("OPT","objective_function")
        except:
            print "ERROR in optimization_algorithm.py: You must indicate the objective function."
            opt.errorexit()
        if self._objective_function not in ["QM_MM_Loss","Physical_Properties_Loss","PhysProp_QMMM_Loss"]:
            print "ERROR in optimization_algorithm.py: Your indicated objective function %s is not in \"QM_MM_Loss\", \"Physical_Properties_Loss\" or \"PhysProp_QMMM_Loss\"."
            opt.errorexit()

    def get_optimization_problem(self):
        """ return the optimization problem """
        return self._optimization_problem

    def get_step_length_control(self):
        """ return the step length control """
        return self._step_length_control

    def get_number_of_iterations(self):
        """ return the number of iterations """
        return self.__number_of_iterations

    def get_tolerance(self):
        """ return the tolerance """
        return self.__tolerance

    def set_tolerance(self, tolerance):
        """ set the tolerance """
        self.__tolerance = tolerance

    def get_number_of_iteration(self):
        """ return the number of the actual iteration """
        return self._number_of_iteration

    def get_iteration(self):
        """ return the actual iteration """
        return self._x
 
    def set_iteration(self, x):
        """ sets the iteration """
        self._x = x

    def get_dimension(self):
        """ return the dimension """
        return self._optimization_problem.get_dimension()

    def execute_optimization_step(self, last_step):
        """ overloaded method to perform optimization step with specific algorithm """
	pass

    def execute_optimization_algorithm(self):
        """ executes the optimization algorithm """

        last_step = False

        for i in range(self.__number_of_iterations):

            self._number_of_iteration = i
   
            # optimum reached?
            if self._optimization_problem.is_optimal(self._x, self.__tolerance):
               last_step = True

            # execute optimization step
            succ = self.execute_optimization_step(last_step)

            # optimization step successful?
            if not succ:
               print "Optimization step was not successful."
               opt.errorexit()

            # finish, if only function evaluation was desired
            if self.__workflow == "feval":
               print "Function evaluation completed."
               break

    def __del__(self):
          """
          destructor
          """	
	  del self	
