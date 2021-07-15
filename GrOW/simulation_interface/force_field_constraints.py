#     class Force_Field_Constraints, 2014/05/19 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		force_field_constraints.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		19-05-2014
#
#     Description:	class Force_Field_Constraints
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
# 			opt_class.py
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

from generic_optimization_problem.box_constraints import Box_Constraints

import math
import opt_class
import string

opt = opt_class.Optimization()

# ***********************************************************************
# class Force_Field_Constraints
# ***********************************************************************

class Force_Field_Constraints(Box_Constraints):
       
    def __init__(self, dimension, initial_guess, config):
          """ constructor """

          Box_Constraints.__init__(self, dimension, initial_guess)

	  self.__config = config
          self.__multicharge = False
          
	  try:
        	self.__param_number = self.__config.get("OPT", "param_number")
 	  except:
   		print "ERROR in force_field_constraints: You must indicate the number of each kind of parameter indicated in your parameter file."
   		opt.errorexit()

 	  self.__param_number = string.split(self.__param_number)
       
          for i in range(len(self.__param_number)):
              self.__param_number[i] = int(self.__param_number[i]) 
  
 	  try:
   		self.__param_order = self.__config.get("OPT", "param_order")
 	  except:
   		print "ERROR in force_field_constraints: You must indicate the order of the parameters in your parameter file."
   		opt.errorexit()

 	  self.__param_order = string.split(self.__param_order)

 	  if len(self.__param_number) != len(self.__param_order):
   	     print "ERROR in force_field_constraints: The length of 'param_number' (%d) does not coincide with the length of 'param_order' (%d)." % (len(self.__param_number), len(self.__param_order))
   	     opt.errorexit()

	     charge_test_list = []
             for i in range(dimension):
                 charge_test_list.append("charge_%d" % (i + 1))

             for param in self.__param_order:
                 if param != "epsilon" and param != "sigma" and param != "dipole" and param != "quadrupole" and param != "charge" and param not in charge_test_list:
                    print "%s is not a valuable force field parameter to optimize. Please choose between 'epsilon', 'sigma', 'dipole', 'quadrupole', and 'charge' (or 'charge_i')!" % (param)
                    opt.errorexit()
                 if param in charge_test_list:
                    self.__multicharge = True
                    
          self.__parind_list = []
          
          ind = 1
          for i in range(len(self.__param_number)):
              ind_list = []
              number = self.__param_number[i]
              for k in range(number):
                  ind_list.append(ind)
                  ind += 1
              self.__parind_list.append(ind_list)
          
          self.__init_boundary(initial_guess)

    def get_param_number(self):
        """ returns the vector of the numbers of parameters """
        return self.__param_numbers

    def get_param_order(self):
        """ returns the vector of the order (names) of parameters """
        return self.__param_order
    
    def get_parind_list(self):
        """ returns the parameter index list """
        return self.__parind_list
 
    def get_multicharge(self):
        """ returns whether the system contains ions or not """
	return self.__multicharge

    def __admissible_domain(self, x, boundary):
       """ help function: computes max and min vector for force field parameters (admissible domain) """

       x_div = []

       k = 0
       for i in range(len(self.__param_order)):
           x_div.append([])
           for j in range(self.__param_number[i]):
               x_div[i].append(x[k])
               k += 1      

       if "mu_L" in self.__param_order:
          nof_params = 2
       else:
          nof_params = len(self.__param_order)

       for i in range(nof_params):
           for j in range(self.__param_number[i]):
               if x_div[i][j] >= 0:
                  self._max_vector.append(x_div[i][j] * (1 + boundary[i] / 100))
                  self._min_vector.append(x_div[i][j] * (1 - boundary[i] / 100))
               else:
                  self._max_vector.append(x_div[i][j] * (1 - boundary[i] / 100))
                  self._min_vector.append(x_div[i][j] * (1 + boundary[i] / 100))

       if "mu_L" in self.__param_order:
          # boundaries for mu (chemical potential)
	  self._min_vector.append(0.0)
          avogadro = 6.02214179
          debye_coulomb = 3.33564
          epsilon_0 = 8.8542
          mu_red_factor = pow(debye_coulomb, 2) / (4 * math.pi * epsilon_0) * avogadro * pow(10, -1)
          try:
              multipole = self.__config.get("OPT","multipole")
          except:
              print "ERROR in force_field_constraints: You must indicate the multipole type."
              opt.errorexit()
          if multipole not in ["dipole","quadrupole"]:
             print "ERROR in force_field_constraints: The variable 'multipole' must be either 'dipole' or 'quadrupole'. %s is not a valid value." %(multipole)
             opt.errorexit()
          if multipole == "quadrupole":
             self._max_vector.append(4 * self._min_vector[0] * pow(self._min_vector[1], 5) / mu_red_factor)
          else:
             self._max_vector.append(4 * self._min_vector[0] * pow(self._min_vector[1], 3) / mu_red_factor)
          # boundaries for L (elongation)
          self._min_vector.append(0.0)
          L_red_factor = 1.0
          self._max_vector.append(0.8 * self._min_vector[1] / L_red_factor)         

    def __init_boundary(self, x):
	""" sets lower and upper bound for loss function's feasible set """

   	try:
     	 	boundary = self.__config.get("OPT", "boundary")
   	except:
      		print "ERROR in force_field_constraints: You must indicate the boundary conditions for your force field parameters."
      		opt.errorexit()

   	boundary = string.split(boundary)
   	boundary = opt.float_vector(boundary, check=True)
   	if type(boundary) != list:
      	   print "ERROR in force_field_constraints: The indicated boundary '%s' is not numeric." % (boundary)
           opt.errorexit()

   	if len(boundary) != len(self.__param_order):
      	   print "ERROR in force_field_constraints: The length of your boundary vector (%d) does not correspond to the length of your parameter order vector (%d)." % (len(boundary), len(self.__param_order))
      	   opt.errorexit()

  	opt.setenv("BOUNDARY", repr(boundary))

        self.__admissible_domain(x, boundary)

    def get_rest_charges(self, x): 
        """ finds out the rest charges summing up to the correct total charges (constraint sum_i q_i = Q) """

        try:
		total_charge = self.__config.get("OPT", "total_charge")
    	except:
		print "ERROR in force_field_constraints: You must indicate the total charge of the molecule simulated."
       		opt.errorexit()

    	total_charge = string.split(total_charge)

    	for i in range(len(total_charge)):
            try:
	   	total_charge[i] = float(total_charge[i])
            except:
	  	print "ERROR in force_field_constraints: Your indicated charge in 'total_charge' (%s) is not a numeric value." % (total_charge[i])
           	opt.errorexit()

        try:
	    stoichiometry = self.__config.get("OPT","stoichiometry")
        except:
	    print "ERROR in force_field_constraints: You must indicate the stoichiometry of the molecule simulated."
            opt.errorexit()

        stoichiometry = stoichiometry.split(",")
    
        if len(stoichiometry) != len(total_charge):
	   print "The length of total_charge is not equal to the length of stoichiometry."
	   opt.errorexit()
    
        for i in range(len(stoichiometry)):
	    stoichiometry[i] = string.split(stoichiometry[i])
  
        for i in range(len(stoichiometry)):
             for j in range(len(stoichiometry[i])):

		try:
	   		stoichiometry[i][j] = float(stoichiometry[i][j])
        	except:
	  		print "Your indicated value in 'stoichiometry' (%s) is not a numeric value." %(stoichiometry[i][j])
          		opt.errorexit()

    	charge_name_list = []
    	for elem in self.__param_order:
            if elem.startswith("charge"):
               charge_name_list.append(elem)

    	if len(total_charge) != len(charge_name_list):
       	   print "ERROR in force_field_constraints: The length of 'total_charge' (%d) is not equal to the number of charges indicated in 'param_order' (%d)." % (len(total_charge), len(charge_name_list))
       	   opt.errorexit()
 
        charge_test_list = []
        for i in range(self._dimension):
          charge_test_list.append("charge_%d" %(i+1))

        list_of_charge_lists = []
        list_of_index_lists = []
        for i in range(len(charge_name_list)):
            if self.__multicharge == True:
               charge_name = "charge_%d" % (i + 1)
            else:
               charge_name = "charge"

        lists = list(opt.get_all_kind_of(charge_name, x, self.__param_order, self.__parind_list))
        list_of_charge_lists.append(lists[0])
        list_of_index_lists.append(lists[1])

        opt.setenv("CHARGES", repr(list_of_charge_lists))
        opt.setenv("CHARGEINDS", repr(list_of_index_lists))
        opt.setenv("TOTALCHARGE", repr(total_charge))
        opt.setenv("STOICHIOMETRY", repr(stoichiometry))

        list_of_charge_lists = []
        
        for index_list in list_of_index_lists:
            charge_list = []
            for index in index_list:
                charge_list.append(x[index])   
            list_of_charge_lists.append(charge_list)

        for i in range(len(stoichiometry)):
 	    if len(stoichiometry[i])-1 != len(list_of_charge_lists[i]):
	       print "The length of 'stoichiometry' (%s) is not equal to the number of charges indicated in 'param_order'." %(i+1)
	       opt.errorexit()
       
        rest_charge = []
         
        for i in range(len(list_of_charge_lists)):

	    stoichiometry_sum =[]
	
	    z = len(stoichiometry[i])
	    for j in range(len(list_of_charge_lists[i])):
	        n = stoichiometry[i][j]*list_of_charge_lists[i][j]
	        stoichiometry_sum.append(n)
	    
	    rest_charge_c = (total_charge[i]-sum(stoichiometry_sum))/stoichiometry[i][z-1]
            rest_charge.append(rest_charge_c)

        return rest_charge

    def __del__(self):
          """
          destructor
          """	
	  del self	
