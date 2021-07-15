#     class Steepest_Descent, 2014/05/14 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		steepest_descent.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		14-05-2014
#
#     Description:	class Steepest_Descent
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
# 			
#     Output:	        
#
#     Imported:		armijo_step_length_control.py
# 			opt_class.py
#                       optimization_algorithm.py
# 			optimization_problem.py
# 			
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

import armijo_step_length_control
from optimization_algorithm import Optimization_Algorithm

import opt_class
opt = opt_class.Optimization()

import os
import sys

# ***********************************************************************
# class Steepest_Descent
# ***********************************************************************

class Steepest_Descent(Optimization_Algorithm):
       
    def __init__(self, optimization_problem, armijo_step_length_control, number_of_iterations, config):
        """ constructor """
          
        self.__armijo = armijo_step_length_control

        Optimization_Algorithm.__init__(self, optimization_problem, self.__armijo, number_of_iterations, config)

    def execute_optimization_step(self, last_step):
        """ perform optimization step with steepest descent algorithm """

        # iteration
        F_new = None
        x = list(self.get_iteration())
        loop = self.get_number_of_iteration()
        opt.setenv("loop", repr(loop + 1))

        # gradient
        opt.setenv("name", "g.%d" % (loop))
        F = self._optimization_problem.get_objective_function()
        opt.setenv("DIM", repr(F.get_dimension()))
        gradient = list(F.get_gradient(x))

        # optimization finished?
        if last_step:
            print "Optimization finished."
            sys.exit()

        if opt.norm(gradient) <= 1.0:
            opt.setenv("DERIV_NORM","no")
    
        # descent = negative gradient
        descent = []
        for i in range(len(gradient)):
            descent.append(-gradient[i])

        # Armijo step length
        if F_new == None:
            F_new = F.get_function_value(x)
            print "F = ", F_new

        opt.setenv("name", "a.%d" % (loop))
        [step_length, F_new] = self.__armijo.get_step_length(eval(os.getenv("MINVECTOR")), eval(os.getenv("MAXVECTOR")), x, F_new, gradient, descent)
        print "F = ", F_new
        l = self.__armijo.get_number_of_iteration()

    	# set new iteration
        x_new = []
        for i in range(len(x)):
            if opt.norm(gradient) > 1.0:
                x_new.append(x[i] + step_length * descent[i] / opt.norm(descent))
            else:
                x_new.append(x[i] + step_length * descent[i])
        self.set_iteration(x_new)

        if self._objective_function == "QM_MM_Loss" or self._objective_function == "Physical_Properties_Loss":
            new_dir = os.path.join(os.getenv("GROW_HOME"), "g.%d.0" % (loop + 1))
            try:
                os.stat(new_dir)
            except:
                os.makedirs(new_dir)

            os.system("cp -r %s/* %s" % (os.path.join(os.getenv("GROW_HOME"), os.getenv("name") + ".%d" % (l)), new_dir))
        if self._objective_function == "PhysProp_QMMM_Loss":
            new_dir_QM_MM = os.path.join(os.getenv("GROW_HOME_QM_MM"), "g.%d.0" % (loop + 1))
            try:
                os.stat(new_dir_QM_MM)
            except:
                os.makedirs(new_dir_QM_MM)

            os.system("cp -r %s/* %s" % (os.path.join(os.getenv("GROW_HOME_QM_MM"), os.getenv("name") + ".%d" % (l)), new_dir_QM_MM))

            new_dir_Physical_Properties = os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), "g.%d.0" % (loop + 1))
            try:
                os.stat(new_dir_Physical_Properties)
            except:
                os.makedirs(new_dir_Physical_Properties)

            os.system("cp -r %s/* %s" % (os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), os.getenv("name") + ".%d" % (l)), new_dir_Physical_Properties))

        print "New iteration: ", x_new
        return True

    def __del__(self):
        """
        destructor
        """	
        del self	
