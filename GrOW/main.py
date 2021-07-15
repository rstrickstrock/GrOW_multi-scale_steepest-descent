#     main GROW control script, 2014/06/06 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		main.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		06-06-2014
#
#     Description:	main GROW control script
#
#     Usage:		main.py <configuration_file>
#
#     Arguments:	configuration_file: name of configuration file
#
#     Options:	       
# 			
#     Output:	        
#
#     Imported:		opt_class.py
#
#     Called:		all required subscripts
#
#     Modifications:
#
# 	$Log: md_opt$
#
#     end_doc

# test


# *****************************************************************************
# Title
# *****************************************************************************

print ""
print " ================================================== "
print " || GROW V. 2.0 - The modern way of optimization || "
print " ================================================== "
print ""
print " Author/Concept:        Marco Huelsmann "
print " Additional Authors:    Andreas Kraemer"
print "                        Astrid Maass"
print "                        Janina Hemmersbach"
print "                        Doron Dominic Heinrich"
print "                        Markus Huber"
print "                        Karl N. Kirschner"
print "                        Sonja Kopp"
print "                        Ottmar Kraemer-Fuhrmann"
print "                        Thomas J. Mueller"
print "                        Thorsten Koeddermann"
print "                        Robin Strickstrock"
print "                        Andre Tissen"
print ""
print " When using GROW, please reference: "
print "   M. Huelsmann, T. Koeddermann, J. Vrabec, D. Reith "
print "   GROW: A Gradient-based Optimization Workflow for the Automated "
print "   Development of Molecular Models"
print "   Computer Physics Communications 181 (2010), 499-513 "
print ""

# ***********************************************************************
# required python modules
# ***********************************************************************

from generic_optimizer.armijo_step_length_control import Armijo_Step_Length_Control
#from generic_optimizer.conjugate_gradient import Conjugate_Gradient
#from generic_optimizer.spagrow import Spagrow
from generic_optimizer.steepest_descent import Steepest_Descent
#from generic_optimizer.cosmos import CoSMoS
from simulation_interface.force_field_constraints import Force_Field_Constraints
from simulation_interface.molecular_simulation_optimization_problem import Molecular_Simulation_Optimization_Problem
from simulation_interface.physical_properties_loss_new import Physical_Properties_Loss 
#from simulation_interface.physical_properties_loss import Physical_Properties_Loss
#from simulation_interface.qm_mm_loss import QM_MM_Loss
#from simulation_interface.spectrum_loss import Spectrum_Loss

from utilities.io import IO
i_o = IO()

from utilities.system import System
sy = System()

from utilities.trace import Trace
tr = Trace()

import os
import os.path
import string
import sys
import time


# *****************************************************************************
# System check function
# *****************************************************************************
def check_system(objective_function):
   """ checks general system specific variables and paths """

   GROW_HOME = os.getenv("GROW_HOME")
   #module = os.getenv("MODULE")
   #print "MODULE: %s" %(module)
   #exit()

   if GROW_HOME == None:
      print "You must set the global environment variable 'GROW_HOME'!" 
      tr.errorexit()
       
   try:
      os.stat(GROW_HOME)
   except:
      print "Your indicated 'GROW_HOME' directory (%s) does not exist." %(GROW_HOME)
      tr.errorexit()

   if objective_function == "PhysProp_QMMM_Loss":
      GROW_HOME_QM_MM = os.path.join(GROW_HOME,"QM_MM")
      GROW_HOME_PHYSICAL_PROPERTIES = os.path.join(GROW_HOME,"Physical_Properties")

      sy.setenv("GROW_HOME_QM_MM",GROW_HOME_QM_MM)
      sy.setenv("GROW_HOME_PHYSICAL_PROPERTIES",GROW_HOME_PHYSICAL_PROPERTIES)

      try:
         os.stat(GROW_HOME_QM_MM)
      except:
         os.mkdir(GROW_HOME_QM_MM)
         
      #else:
      #   print "The QM_MM home directory %s already exists. To prevent data loss the optimization will exit." %(GROW_HOME_QM_MM)
      #   tr.errorexit()

      try:
         os.stat(GROW_HOME_PHYSICAL_PROPERTIES)
      except:
         os.mkdir(GROW_HOME_PHYSICAL_PROPERTIES)

      #else:
      #   print "The PHYSICAL_PROPERTIES home directory %s already exists. To prevent data loss the optimization will exit." %(GROW_HOME_PHYSICAL_PROPERTIES)
      #   tr.errorexit()

    

# ***********************************************************************
# reading the command line
# ***********************************************************************

# # Arguments of system call
args = sys.argv

# initialize standard variables
showhelp = 'n'                # default: no help to be shown
showinfo = 'n'                # default: no info to be shown
wizard = 'y'                  # default: no wizard to complete input on
usedefaults = 'y'             # default: use 
verbose = 7                   # default: middle verbosity
config_file = "not_defined"   

counter = 1
while counter < len(args):
   this_argument = args[counter]
   if this_argument[0] == '-':
      # parsing options
      done = 'n'
      if this_argument == "-h":
         showhelp = 'y'
         done = 'y'
      if this_argument == "-i":
         showinfo = 'y'
         done = 'y'
      if this_argument == "-nw":
         wizard = 'n'
         done = 'y'
      if this_argument == "-nd":
         usedefaults = 'n'
         done = 'y'
      if this_argument == "-v":
         if counter + 1 < len(args):
            try:
            	next_argument = int(args[counter+1])
            except:
		print "Warning: option -v was not followed by an integer."
                print "         verbose level will remain at 7 (default)."
                print ""
                done = 'y'
            else:
                verbose = next_argument
            counter = counter + 1
            done = 'y'

	 else:
            print "Warning: option -v was not followed by an integer."
            print "         verbose level will remain at 7 (default)."
            print ""
            done = 'y'
      if done == 'n':
         print "Warning: option %s is not known and will be" %(this_argument)
         print "         neglected."
         print ""
      counter = counter + 1  # go to the next argument in the command line
   else:
      # testing for config file
      config_file = this_argument
      try:
         os.stat(config_file)
      except:
         print "Warning: The following non-option-input in the command line"
         print "         was not the expected configuration file and will be"
         print "         neglected: %s" %(config_file)
         print ""
         config_file = "not_defined" 
      counter = counter + 1  # got to the next argument in the command line

sy.setenv("VERBOSITY",repr(verbose))


# ***********************************************************************
# reading the input
# ***********************************************************************
try:
	os.stat(args[1])
except OSError: 
	print "ERROR in main.py: Your indicated configuration file '%s' does not exist." % (args[1])
        tr.errorexit()
except:
	print "ERROR in main.py: You must indicated the configuration file."
        tr.errorexit()

# # Reading <configuration_file>
config_file = os.path.abspath(args[1])
config = i_o.read_config_file(config_file)
sy.setenv("CONFIGFILE", config_file)

try:
	objective_function = config.get("OPT", "objective_function")
except:
	print "ERROR in main.py: You must indicate the objective function."
	tr.errorexit()

if objective_function not in ["Physical_Properties_Loss", "QM_MM_Loss", "Spectrum_Loss","PhysProp_QMMM_Loss"]:
	print "ERROR in main.py: The only valid objective functions are Physical_Properties_Loss and QM_MM_Loss and Spectrum_Loss and PhysProp_QMMM_Loss . '%s' is not valid." % (objective_function)
	tr.errorexit()

# Reading the config file
# -----------------------

# ***********************************************************************
# checking the global system
# ***********************************************************************

check_system(objective_function)


# ***********************************************************************
# defining the optimization problem and algorithm
# ***********************************************************************


# Initial guess
# -------------


if objective_function == "Physical_Properties_Loss" or objective_function == "QM_MM_Loss":
   try: 
      par = config.get("OPT", "parameters")
   except:
      print "ERROR in main.py: You must indicate a parameter file."
      tr.errorexit()

   try:
      os.stat(par)
   except:
      print "ERROR in main.py: Parameter file %s does not exist." % (par)
      tr.errorexit()
	
   initial_guess = i_o.read_last_parameter(par)

   # Dimension
   # ---------
   dimension = len(initial_guess)


   # objective function

   weights = []
   targets = []

   #if objective_function == "Physical_Properties_Loss":
   #   obj_fct = Physical_Properties_Loss(dimension, weights, targets, config)
#   if objective_function == "QM_MM_Loss":
#      obj_fct = QM_MM_Loss(dimension, weights, targets, config)
   #elif objective_function == "Spectrum_Loss":
   #   obj_fct = Spectrum_Loss(dimension, weights, targets, config)

   ### changed: only one class with both functionalities of physical_properties_loss and qm_mm_loss called physical_properties_loss_new
   ### in config the objective_function in [OPT] is now used within the class and not to decide which class should be used
   ### -> NO USABILITY CHANGES!

   obj_fct = Physical_Properties_Loss(dimension, weights, targets, config)

#   sy.setenv("WEIGHTS", repr(obj_fct.get_weights()))


   # force field constraints

   try:
      constraints = config.get("OPT", "constraints")
   except:
      print "ERROR in main.py: You must indicate the constraints."
      tr.errorexit()

   if constraints != "Force_Field_Constraints":
      print "ERROR in main.py: The only valid constraints are Force_Field_Constraints. '%s' is not valid." % (constraints)
      tr.errorexit()

   ff_constr = Force_Field_Constraints(dimension, initial_guess, config)
   [min_vector, max_vector] = ff_constr.get_boundary()
   sy.setenv("MINVECTOR", repr(min_vector))
   sy.setenv("MAXVECTOR", repr(max_vector))

   # optimization problem
   try:
      optimization_problem = config.get("OPT", "optimization_problem")
   except:
      print "ERROR in main.py: You must indicate the optimization problem."
      tr.errorexit()

   if optimization_problem != "Molecular_Simulation_Optimization_Problem":
      print "ERROR in main.py: The only valid optimization problem Molecular_Simulation_Optimization_Problem. '%s' is not valid." % (optimization_problem)
      tr.errorexit()

   opt_prob = Molecular_Simulation_Optimization_Problem(dimension, obj_fct, ff_constr)

   # optimization algorithm
   try:
      method = config.get("OPT", "method")
   except:
      print "ERROR in main.py: You must indicate the optimization method you want to use."
      tr.errorexit()
         
   if method == "Steepest_Descent" or method == "Conjugate_Gradient":
      step_length_obj = Armijo_Step_Length_Control(opt_prob, config)
      number_of_iterations = 100
      if method == "Steepest_Descent":
         opt_algo = Steepest_Descent(opt_prob, step_length_obj, number_of_iterations, config)
      else:    
         opt_algo = Conjugate_Gradient(opt_prob, step_length_obj, number_of_iterations, config)

      opt_algo.execute_optimization_algorithm()

   if method == "cosmos":    
      step_length_obj = None
      number_of_iterations = None
      opt_algo = CoSMoS(opt_prob, step_length_obj, number_of_iterations, config)
      opt_algo.execute_optimization_algorithm()

   if method == "Spagrow":    
      number_of_iterations = 100
      opt_algo = Spagrow(opt_prob, number_of_iterations, config)
      opt_algo.execute_optimization_algorithm()



if objective_function == "PhysProp_QMMM_Loss":
   try: 
      par = config.get("OPT", "parameters")
   except:
      print "ERROR in main.py: You must indicate a parameter file."
      tr.errorexit()

   try:
      os.stat(par)
   except:
      print "ERROR in main.py: Parameter file %s does not exist." % (par)
      tr.errorexit()
	
   initial_guess = i_o.read_last_parameter(par)

   # Dimension
   # ---------
   dimension = len(initial_guess)

   # objective function
   weights_QM_MM = []
   weights_Physical_Properties = []
   weights = [weights_QM_MM, weights_Physical_Properties]

   targets_QM_MM = []
   targets_Physical_Properties = []
   targets = [targets_QM_MM, targets_Physical_Properties]

   obj_fct = Physical_Properties_Loss(dimension, weights, targets, config)

   sy.setenv("WEIGHTS_QM_MM", repr(obj_fct.get_weights_QM_MM()))
   sy.setenv("WEIGHTS_PHYSICAL_PROPERTIES", repr(obj_fct.get_weights_Physical_Properties()))


   # force field constraints
   try:
      constraints = config.get("OPT", "constraints")
   except:
      print "ERROR in main.py: You must indicate the constraints."
      tr.errorexit()

   if constraints != "Force_Field_Constraints":
      print "ERROR in main.py: The only valid constraints are Force_Field_Constraints. '%s' is not valid." % (constraints)
      tr.errorexit()

   ff_constr = Force_Field_Constraints(dimension, initial_guess, config)
   [min_vector, max_vector] = ff_constr.get_boundary()
   sy.setenv("MINVECTOR", repr(min_vector))
   sy.setenv("MAXVECTOR", repr(max_vector))

   # optimization problem

   try:
      optimization_problem = config.get("OPT", "optimization_problem")
   except:
      print "ERROR in main.py: You must indicate the optimization problem."
      tr.errorexit()

   if optimization_problem != "Molecular_Simulation_Optimization_Problem":
      print "ERROR in main.py: The only valid optimization problem Molecular_Simulation_Optimization_Problem. '%s' is not valid." % (optimization_problem)
      tr.errorexit()

   opt_prob = Molecular_Simulation_Optimization_Problem(dimension, obj_fct, ff_constr)

   # optimization algorithm
   try:
      method = config.get("OPT", "method")
   except:
      print "ERROR in main.py: You must indicate the optimization method you want to use."
      tr.errorexit()
         
   if method == "Steepest_Descent" or method == "Conjugate_Gradient":
      step_length_obj = Armijo_Step_Length_Control(opt_prob, config)
      number_of_iterations = 100
      if method == "Steepest_Descent":
         opt_algo = Steepest_Descent(opt_prob, step_length_obj, number_of_iterations, config)
      else:
         print "method \"Conjugate Gradient\" is not tested (or implemented?), yet. Change to  \"Steepes Descent\" (or check on your own in main.py) ."
         tr.errorexit()    
         opt_algo = Conjugate_Gradient(opt_prob, step_length_obj, number_of_iterations, config)

      opt_algo.execute_optimization_algorithm()

   if method == "cosmos":    
      step_length_obj = None
      number_of_iterations = None
      opt_algo = CoSMoS(opt_prob, step_length_obj, number_of_iterations, config)
      opt_algo.execute_optimization_algorithm()

   if method == "Spagrow":    
      number_of_iterations = 100
      opt_algo = Spagrow(opt_prob, number_of_iterations, config)
      opt_algo.execute_optimization_algorithm()



