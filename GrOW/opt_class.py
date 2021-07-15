import os
import os.path
import sys
import traceback
import string
import math
import random as rndm
import time
import efficiency

class Optimization:
       
    def __init__(self):
          """ constructor """
          self.obj = None

          if os.getenv("efficiency") != None:
             self.par_red_list = []
             self.actual_parameters = []
             self.basis = []
             self.eps_ref = None
             self.sigma_ref = None
             self.epsilon_ref_tilde_list = []
             self.sigma_ref_tilde_list = []

    def errorexit(self,tool=False):
       """gives out general error message in case of severe error"""
       if tool == False:
          print "Severe error, aborting."
       else:
          print "Simulation tool error, aborting" 
       try: 
	 	os.stat(os.path.join(os.getenv("TMPPATH"),"gradient.txt"))
       except:
	        pass
       else:
       		os.system("cp %s %s"%(os.path.join(os.getenv("TMPPATH"),"gradient.txt"),os.path.join(os.getenv("OPTOUTPATH"),"gradient.txt")))
       if tool == False:
          self.delete()

       at = time.localtime()
       dt = os.times()

       verbose = int(os.getenv("VERBOSITY"))

       if verbose > 5:
          print "Break Time: %s.%s.%s, %s:%s:%s" %(at[2],at[1],at[0],at[3],at[4],at[5])
       if verbose > 9:
          print "User time:        ", dt[0]
          print "System time       ", dt[1]
          print "Kind user time:   ", dt[2]
          print "Kind system time: ", dt[3]
          print "Elapse real time: ", dt[4]

       sys.exit(1)

    ### Methods    

    def setenv(self,env,value):
       """
       sets the environment variable env to value value
       """

       if value != None:
          os.environ[env] = value
       else:
          if env in os.environ:
             del os.environ[env]

       return


    def read_config_file(self,config_file):
       """
       reads configuration file "config_file" and returns a config object
       """

       try:
	  import ConfigParser
       except:
	  print "No module named ConfigParser but required."
	  print "Abort."
	  self.errorexit()

       config = ConfigParser.ConfigParser()

       f_ptr = open(config_file)

       config.readfp(f_ptr)
       f_ptr.close()
 
       return config

    def check_source_code(self,source_file_list):
       """
       checks the existence of the source code files in current dir and execute dir
       returns updated source_file_list with absolute paths
       """
        
       for i in range(len(source_file_list)) :
          if source_file_list[i] != []:
           found = False
           try:
	      os.stat(os.path.join(os.getcwd(),source_file_list[i]))
           except:
	      try:
	   	   os.stat(os.path.join(os.getenv("EXEDIR"),source_file_list[i]))
              except:
                   try:
			os.stat(os.path.join(os.getenv("BINDIR"),source_file_list[i]))
		   except:
                        pythonpath = string.split(os.getenv("PYTHONPATH"),":")
                        for path in pythonpath:
                            try:
				os.stat(os.path.join(path,source_file_list[i]))
                            except:
	    	        	pass
                            else:
				source_file_list[i] = os.path.join(path,source_file_list[i])
                                found = True
                   else:
			source_file_list[i] = os.path.join(os.getenv("BINDIR"),source_file_list[i])
                        found = True
              else:
                source_file_list[i] = os.path.join(os.getenv("EXEDIR"),source_file_list[i])
                found = True
           else:
              source_file_list[i] = os.path.join(os.getcwd(),source_file_list[i])
              found = True
           if not found:
              print "Can't find program %s, but I need it." %(source_file_list[i])
	      self.errorexit()

       return source_file_list

    def admissible_domain(self,x,boundary,param_order,param_number,dim):
       """ computes max and min vector for force field parameters (admissible domain) """


       max_vector = []
       min_vector = []

       x_div = []

       k = 0
       for i in range(len(param_order)):
           x_div.append([])
           for j in range(param_number[i]):
               x_div[i].append(x[k])
               k += 1      

       for i in range(len(param_order)):
           for j in range(param_number[i]):
               if x_div[i][j] >= 0:
                  max_vector.append(x_div[i][j]*(1+boundary[i]/100))
                  min_vector.append(x_div[i][j]*(1-boundary[i]/100))
               else:
                  max_vector.append(x_div[i][j]*(1-boundary[i]/100))
                  min_vector.append(x_div[i][j]*(1+boundary[i]/100))

       return [max_vector,min_vector]

    def read_dim(self):
       """ reads the dimension from the parameter file """

       param_file = os.getenv("PARAMFILE")

       x = os.popen("tail -1 %s" %(param_file)).readlines()
       x = string.split(x[0]," ")

       return len(x)

    def reduce_korr_parameters(self,x):
        """ produces L*, mu* """

        epsilon = x[0]
        sigma = x[1]
        mu = x[2]
        L = x[3]


        mu_red_factor = float(os.getenv("mu_red_factor"))
        if os.getenv("MULTIPOLE") == "quadrupole":
           mu_red = mu/(epsilon*pow(sigma,5))*mu_red_factor
        else:
           mu_red = mu/(epsilon*pow(sigma,3))*mu_red_factor

        L_red_factor = float(os.getenv("L_red_factor"))
        L_red = L/sigma*L_red_factor

        return [mu_red,L_red]

    def introduce_uncertainties(self,rho,s_rho,vapor,s_vap):
        """ introduces artificial uncertainties into results of correlation function """


        s_rho_percent = rndm.uniform(0.0,s_rho)
        s_vap_percent = rndm.uniform(0.0,s_vap)

        s_rho = rndm.choice([-1,1])*s_rho_percent*rho
        s_vap = rndm.choice([-1,1])*s_vap_percent*vapor
        
        rho_biased = rho + s_rho
        vapor_biased = vapor + s_vap

        f = open(os.path.join(os.getenv("OUTPATH"),"uncertainties.txt"),"w")
        f.write(repr(s_vap)+"\n")
        f.write(repr(s_rho)+"\n")
        f.close()

        return [rho_biased,vapor_biased]

    def efficiency_parameters(self,name):
        """ sets class variables in the case of efficiency """      

        loop = int(os.getenv("loop"))


        if name.startswith("original"):
           eff_obj = efficiency.Efficiency()
           self.actual_parameters = self.read_last_parameter(os.getenv("PARAMFILE"))
           dim = len(self.actual_parameters)
           [self.par_red_list,self.epsilon_ref,self.sigma_ref] = eff_obj.reduce_parameters(self.actual_parameters,dim)
           self.basis = []
        if name.startswith("gradient"):
           y = self.read_last_parameter(os.getenv("PARAMFILE"))
           self.basis.append(self.vector_sum(y,self.mult_vector(-1.0,self.actual_parameters)))

           if len(self.basis) == len(self.actual_parameters):
              self.test_independence(self.basis,loop-1)

        return

    def get_rest_charges(self,x):
        """ finds out the rest charges summing up to the correct total charges """

        list_of_index_lists = eval(os.getenv("CHARGEINDS"))
        total_charge = eval(os.getenv("TOTALCHARGE"))

        list_of_charge_lists = []
        
        for index_list in list_of_index_lists:
            charge_list = []
            for index in index_list:
                charge_list.append(x[index])   
            list_of_charge_lists.append(charge_list)
       
        rest_charge = []
         
        for i in range(len(list_of_charge_lists)):
            rest_charge.append(total_charge[i]-sum(list_of_charge_lists[i]))

        return rest_charge


    def change_units(self, pars, order, parind_list):
        """ converts sigma from nm to A and epsilon from kJ/mol to K """

        for i in range(len(pars)):
          try:
            if self.which_parameter(order,parind_list,i+1) == "sigma":
               pars[i] = pars[i]*10
            if self.which_parameter(order,parind_list,i+1) == "epsilon":
               pars[i] = pars[i]*1000/(6.02214*1.38066)
          except:
            pass

        return pars

    def start_simulation(self,program,md_aux,param_index,loop=None,NMOL=None,name=None):
        """ starts a molecular simulation run """

        tmppath = os.getenv("TMPPATH")

        if program == "korr":
           try:
		os.stat(os.path.join(tmppath,name))
           except:
		os.makedirs(os.path.join(tmppath,name))

           x = self.read_last_parameter(os.getenv("PARAMFILE"))
      
           if os.getenv("TRANGE_ERROR") != None:
              T_range_err = eval(os.getenv("TRANGE_ERROR"))
           else:
              T_range_err = [float(os.getenv("TEMPERATURE"))]

           [mu_red,L_red] = self.reduce_korr_parameters(x)         

           epsilon = x[0]
           sigma = x[1]

	   # Ausgleich moeglicher Rundungsfehler
           if os.getenv("MULTIPOLE") == "quadrupole":
              if mu_red > 12.0:
                 mu_red = 12.0
           else:
              if mu_red > 4.0:
                 mu_red = 4.0
           if mu_red < 0.0:
              mu_red = 0.0
           if L_red > 0.5:
              L_red = 0.5
           if L_red < 0.0:
              L_red = 0.0 

           rho_list = []
           vapor_list = []    

           for T in T_range_err:

               T_red_factor = float(os.getenv("T_red_factor"))
               T_red = T/epsilon*T_red_factor

               R_script = os.getenv("MDPROG")

               cmd = "R --slave < %s --args %f %f %f" %(R_script,mu_red,L_red,T_red)
          
               print "Evaluating correlation function, T = %.1f..." %(T)
               res = os.popen(cmd).readlines()
               for i in range(len(res)):
                   res[i] = float(res[i][4:(len(res[i])-1)])

               T_c = res[0]; rho_c = res[1]; rho_l = res[2]; rho_v = res[3]; p_s = res[4]; w = res[5]; c2 = res[6]; c3 = res[7]

               # density
           
               rho_red_factor = float(os.getenv("rho_red_factor"))
               rho = rho_l/pow(sigma,3)*rho_red_factor

               if os.getenv("WHICH") == "vap":
                  # enthalpy of vaporization (by Clausius-Clapeyron equation

                  if T_red/T_c >= 0.7:
                     help_constant = rho_v
		  else:
                     help_constant = p_s/T_red

                  cc_factor = p_s*T_red*(1/help_constant - 1/rho_l)             
                  ln_ps_derivative = -c2/pow(T_red,2) - 4*c3/pow(T_red,5)
                  vapor_red = cc_factor*ln_ps_derivative
     
                     

                  vapor_red_factor = float(os.getenv("vapor_red_factor"))
                  vapor = vapor_red*epsilon*vapor_red_factor
            
                  

               else:
                  # vapor pressure                 

                  pressure_red_factor = float(os.getenv("pressure_red_factor"))
                  vapor = p_s*epsilon/(pow(sigma,3))*pressure_red_factor

               rho_list.append(rho)
               vapor_list.append(vapor)
         
           s_rho = float(os.getenv("SRHO"))
           s_vap = float(os.getenv("SVAP"))

           if s_rho != 0.0 or s_vap != 0.0:
              for i in range(len(rho_list)):
                  [rho_list[i],vapor_list[i]] = self.introduce_uncertainties(rho_list[i],s_rho,vapor_list[i],s_vap)

           out_file = os.path.join(os.getenv("OUTPATH"),"properties.txt")           

           f = open(out_file,"w")
           
	   for i in range(len(rho_list)):                    
               f.write(repr(vapor_list[i])+"\n")
               f.write(repr(rho_list[i])+"\n")
           f.close()
 

           # copy for consistency reasons (!!! out_file = properties_file !!!)
           os.system("cp %s %s" %(out_file,os.path.join(os.getenv("OUTPATH"),"properties_cp.txt")))

        else:
         eff = os.getenv("efficiency")
       
         if eff == "y":
           if name.startswith("original") or name.startswith("armijo"):
              sim_indic = True
           else:
              sim_indic = False
           self.efficiency_parameters(name)
              


         else:
           sim_indic = True

           

         if program == "yasp" and sim_indic == True:
           try: 
 		os.stat(os.path.join(os.getenv("OUTPATH"),"temp.out.%d" %(int(loop))))
           except:
		pass
	   else:
		os.remove(os.path.join(os.getenv("OUTPATH"),"temp.out.%d" %(int(loop))))
           
           
           self.setenv("sim_prog",md_aux[2])
           par_string = "%d %d %d" %(int(param_index),int(loop),int(NMOL))
           self.setenv("PARSTRING",par_string)
           
         if program != "yasp" and sim_indic == True:
           x = self.read_last_parameter(os.getenv("PARAMFILE"))
           self.setenv("sim_prog",md_aux[0])
           par_string = ""
           for i in range(len(x)):
               if i != len(x)-1:
                  par_string += "%.5f " %(x[i])
               else:
                  par_string += "%.5f" %(x[i])
           if os.getenv("CHARGEOPT") == "True":
              rest_charge_file = os.getenv("PARAMFILE")+"_rest_charges"
              rest_charge = self.read_last_parameter(rest_charge_file)
              for charge in rest_charge:
                  par_string += " %.5f" %(charge)
           self.setenv("PARSTRING",par_string)


         if sim_indic == True:
           self.setenv("name",name)
           dist_script = os.getenv("DISTSCRIPT")
           
           fits = os.getenv("fits")
           properties = eval(os.getenv("PROPERTIES"))
           VLE = os.getenv("VLE")

           if VLE == "y":
              VLE_T_range = eval(os.getenv("VLETRANGE"))
              VLE_T_string = ""
              for i in range(len(VLE_T_range)):
                  if i != len(VLE_T_range)-1:
                     VLE_T_string += repr(VLE_T_range[i])+" "
                  else:
                     VLE_T_string += repr(VLE_T_range[i])
              self.setenv("TVLESTRING",VLE_T_string)
              VLE_P_range = eval(os.getenv("PRANGEVLE"))
	      VLE_P_string = ""
              for i in range(len(VLE_P_range)):
                  if i != len(VLE_P_range)-1:
                     VLE_P_string += repr(VLE_P_range[i])+" "
                  else:
                     VLE_P_string += repr(VLE_P_range[i])
              self.setenv("PVLESTRING",VLE_P_string)

           if os.getenv("transport") == "TRUE":
              T_trans = eval(os.getenv("TTRANS"))
              T_trans_string = ""
              for i in range(len(T_trans)):
                  if i != len(T_trans)-1:
                     T_trans_string += repr(T_trans[i])+" "
                  else:
                     T_trans_string += repr(T_trans[i])
              self.setenv("TTRANSSTRING",T_trans_string)

              P_trans = eval(os.getenv("PTRANS"))
              P_trans_string = ""
              for i in range(len(P_trans)):
                  if i != len(P_trans)-1:
                     P_trans_string += repr(P_trans[i])+" "
                  else:
                     P_trans_string += repr(P_trans[i])
              self.setenv("PTRANSSTRING",P_trans_string)

           if os.getenv("BULKNPT") == "y":
              T_bulk_NPT = eval(os.getenv("NPTTRANGE"))
              T_bulk_NPT_string = ""
              for i in range(len(T_bulk_NPT)):
                  if i != len(T_bulk_NPT)-1:
                     T_bulk_NPT_string += repr(T_bulk_NPT[i])+" "
                  else:
                     T_bulk_NPT_string += repr(T_bulk_NPT[i])
              self.setenv("NPTTSTRING",T_bulk_NPT_string)

              P_bulk_NPT = eval(os.getenv("NPTPRANGE"))
              P_bulk_NPT_string = ""
              for i in range(len(P_bulk_NPT)):
                  if i != len(P_bulk_NPT)-1:
                     P_bulk_NPT_string += repr(P_bulk_NPT[i])+" "
                  else:
                     P_bulk_NPT_string += repr(P_bulk_NPT[i])
              self.setenv("NPTPSTRING",P_bulk_NPT_string)

           if os.getenv("BULKNVT") == "y":
              T_bulk_NVT = eval(os.getenv("NVTTRANGE"))
              T_bulk_NVT_string = ""
              for i in range(len(T_bulk_NVT)):
                  if i != len(T_bulk_NVT)-1:
                     T_bulk_NVT_string += repr(T_bulk_NVT[i])+" "
                  else:
                     T_bulk_NVT_string += repr(T_bulk_NVT[i])
              self.setenv("NVTTSTRING",T_bulk_NVT_string)

              rho_bulk_NVT = eval(os.getenv("NVTRHORANGE"))
              rho_bulk_NVT_string = ""
              for i in range(len(rho_bulk_NVT)):
                  if i != len(rho_bulk_NVT)-1:
                     rho_bulk_NVT_string += repr(rho_bulk_NVT[i])+" "
                  else:
                     rho_bulk_NVT_string += repr(rho_bulk_NVT[i])
              self.setenv("NVTRHOSTRING",rho_bulk_NVT_string)

           if fits == None or (fits == "y" and "density" not in properties) or (fits == "y" and os.getenv("density_fit") == "linear") or VLE == "y":
              if ("vapor" in properties and program == "gromacs") or program == "ms2":
                 self.change_units(par_string,program)
              # call simulation
              res = os.system("%s" %(dist_script))
              if os.getenv("status") or res != 0:
	         self.errorexit(tool=True)
              T_exists = []
              T_total = []
              if os.getenv("VLE") == "y":
                 VLE_T_range = eval(os.getenv("VLETRANGE"))
                 for k in range(len(VLE_T_range)):
                     T_total.append(VLE_T_range[k])
              if os.getenv("transport") == "TRUE":
                  T_trans = eval(os.getenv("TTRANS"))
                  for k in range(len(T_trans)):
                      T_total.append(T_trans[k])
              if os.getenv("BULKNPT") == "y":
                  T_bulk_NPT = eval(os.getenv("NPTTRANGE"))
                  for k in range(len(T_bulk_NPT)):
                      T_total.append(T_bulk_NPT[k])
              if os.getenv("BULKNVT") == "y":
                  T_bulk_NVT = eval(os.getenv("NVTTRANGE"))
                  for k in range(len(T_bulk_NVT)):
                      T_total.append(T_bulk_NVT[k])

              while len(T_exists) != len(T_total):
                 time.sleep(60)
                 for k in range(len(VLE_T_range)):
                     file_to_test = os.path.join(tmppath,repr(T_total[k]),name,"terminatedVLE.txt") 
                     os.stat(file_to_test)
                     try:
			os.stat(file_to_test)
                     except:
			pass
                     else:
                        if file_to_test not in T_exists:
                           T_exists.append(file_to_test)
                 for k in range(len(T_trans)):
                     file_to_test = os.path.join(tmppath,repr(T_total[len(VLE_T_range)+k]),name,"terminatedtrans.txt") 
                     os.stat(file_to_test)
                     try:
			os.stat(file_to_test)
                     except:
			pass
                     else:
                        if file_to_test not in T_exists:
                           T_exists.append(file_to_test)
                 for k in range(len(T_bulk_NPT)):
                     file_to_test = os.path.join(tmppath,repr(T_total[len(VLE_T_range)+len(T_trans)+k]),name,"terminatedNPT.txt") 
                     os.stat(file_to_test)
                     try:
			os.stat(file_to_test)
                     except:
			pass
                     else:
                        if file_to_test not in T_exists:
                           T_exists.append(file_to_test)
                 for k in range(len(T_bulk_NVT)):
                     file_to_test = os.path.join(tmppath,repr(T_total[len(VLE_T_range)+len(T_trans)+len(T_bulk_NPT)+k]),name,"terminatedNVT.txt") 
                     os.stat(file_to_test)
                     try:
			os.stat(file_to_test)
                     except:
			pass
                     else:
                        if file_to_test not in T_exists:
                           T_exists.append(file_to_test)
		
           at = time.localtime()
           dt = os.times()

           verbose = int(os.getenv("VERBOSITY"))

           if verbose > 5:
              print "Actual Time: %s.%s.%s, %s:%s:%s" %(at[2],at[1],at[0],at[3],at[4],at[5])
           if verbose > 9:
              print "User time:        ", dt[0]
              print "System time       ", dt[1]
              print "Kind user time:   ", dt[2]
              print "Kind system time: ", dt[3]
              print "Elapse real time: ", dt[4]

           os.system("rm %s" %(os.path.join(os.getenv("ROOTDIR"),"*")))

        if os.getenv("SUBMISSIONTEST") == "y":
           if res != 0:
              print "Submission test failed."
              self.errorexit()
           else:
              print "Submission test successful."
              sys.exit(1)

        if os.getenv("WORKFLOW") == "simulation":
           print "Simulation performed successfully."
           sys.exit(1)

    def eval_vector(self,vector,check=False):
        """ evaluates the components of a string vector """
 
        for i in range(len(vector)):
            if check == False:
               vector[i] = eval(vector[i])
            else:
               	try:
   			vector[i] = eval(vector[i])
		except:
			return vector[i]
        return vector

    def float_vector(self,vector,check=False):
        """ floats the components of a string vector """
 
        for i in range(len(vector)):
            if check == False:
               vector[i] = float(vector[i])
            else:
               	try:
   			vector[i] = float(vector[i])
		except:
			return vector[i]
        return vector


    def read_last_parameter(self,filename,fl=True):
        """ reads the last parameter from a file """


        x = os.popen("tail -1 %s" %(filename)).readlines()
        x = string.split(x[0]," ")
        x[len(x)-1] = x[len(x)-1][:(len(x[len(x)-1])-1)]
        if fl == True:
           x = self.float_vector(x)

        return x

    def read_last_two_parameters(self,filename,fl=True):
        """ reads the last two parameters from a file """

        dim = int(os.getenv("dim"))

        x = os.popen("tail -2 %s" %(filename)).readlines()
        x1 = string.split(x[0]," ")
        x2 = string.split(x[1]," ")
        x1[len(x1)-1] = x1[len(x1)-1][:(len(x1[len(x1)-1])-1)]
        if fl == True:
           x1 = self.float_vector(x1)
        x2[len(x2)-1] = x2[len(x2)-1][:(len(x2[len(x2)-1])-1)]
        if fl == True:
           x2 = self.float_vector(x2)

        return [x1,x2]

    def read_line(self,filename,listing=False):
        """ reads a file containing one line """
 
        f = open(filename,"r")

        s = f.readlines()[0]
        s = s[:(len(s)-1)]
     
        
        if listing == True:
           s = string.split(s," ")
        f.close() 
        return s

    def get_line_number(self,filename):
        """ returns the number of lines of a file """
 
        f = open(filename,"r")
        s = len(f.readlines())
        f.close() 
        return s

    def read_dict_lines(self,filename,line_number):
        """ reads a certain line of a file """
        f = open(filename,"r")

        s = f.readlines()[line_number]
        s = s[:(len(s)-1)]

        f.close() 
        return s

    def calculate_rho(self,C,T,Tc):
        """ calculates density by Guggenheim fit """
  
        one_third = float(1.0/3.0)

        rho_red = C[0]+C[1]*pow(Tc-T,one_third)+C[2]*(Tc-T)+C[3]*pow(Tc-T,1.5)

        rho_red_factor = float(os.getenv("rho_red_factor"))
        rho = rho_red_factor*rho_red

        return rho

    def calculate_vapor(self,C,T,Tc):
        """ calculates enthalpy of vaporisation by Guggenheim fit """
  
        Tr = T/Tc

        vapor_red = math.exp(C[0])*pow(1-Tr,C[1]+C[2]*Tr+C[3]*pow(Tr,2))

        vapor_red_factor = float(os.getenv("vapor_red_factor"))
        vapor = vapor_red_factor*vapor_red

        return vapor


    def get_out_file(self,program,name=None):
        """ gets the name of the output file of the simulation program """


        if program == "yasp":
           outpath = os.path.join(os.getenv("ORIG_OUTPATH"),os.getenv("ORIG_TEMPERATURE"))
           out_file = os.path.join(outpath,"temp.out.%d" %(int(os.getenv("loop"))))
        elif program == "korr":
           out_file = os.path.join(os.getenv("OUTPATH"),"properties.txt")
        else:
           outpath = os.getenv("ORIG_OUTPATH")
           out_file = os.path.join(outpath,name)
           os.system("touch %s" %(out_file))


        return out_file       
   
           


    def make_properties_file_by_lists(self,outpath,prop):
        """ makes a properties file by lists (e.g. gromacs) """      

        tmppath = os.getenv("TMPPATH")
        aux = eval(os.getenv("auxiliaries"))
        gro_R = aux[0] 

        program = os.getenv("program")

        # density
        if prop == "density" or prop == "sld":
           if prop == "density":
              ROUT = os.path.join(outpath,"Bulk_Density.list")
              sd_file = os.path.join(outpath, "Bulk_Density.sd")
           elif prop == "sld":
              ROUT = os.path.join(outpath,"Density.list")
              sd_file = os.path.join(outpath,"Density.sd")
           self.setenv("ROUT",ROUT)
           self.setenv("AVERAGE_TYPE","LJ")
           res = os.system("R --slave < %s" %(gro_R))
           if res != 0:
              self.errorexit()
           dens_file = os.path.join(tmppath,"average.txt")
           dens_vector = self.read_last_parameter(dens_file)
           dens = self.mean(dens_vector)
           try:
		os.stat(sd_file)
           except:
                drho = self.sd(dens_vector)
           else:
		drho = self.read_last_parameter(sd_file)[0]

           property = dens
           noise = drho     

        # enthalpy of vaporization
        if prop == "vapor":
           ROUT = os.path.join(outpath,"Vapor.list")
           self.setenv("ROUT",ROUT)
           self.setenv("AVERAGE_TYPE","LJ")
           self.setenv("AVERAGE_TYPE","LJ")
           res = os.system("R --slave < %s" %(gro_R))
           if res != 0:
              self.errorexit()
           vapor_file = os.path.join(tmppath,"average.txt")
           vapor_vector = self.read_last_parameter(vapor_file)
           vapor = self.mean(vapor_vector)
           try:
		os.stat(os.path.join(outpath,"Vapor.sd"))
           except:
                dv = self.sd(vapor_vector)
           else:
                sd_file = os.path.join(outpath,"Vapor.sd")
		dv = self.read_last_parameter(sd_file)[0]
       
           property = vapor
           noise = dv

        # (Vapor) Pressure
        if prop == "P" or prop == "pressure":
           if prop == "P":
              ROUT = os.path.join(outpath,"Pressure.list")
           elif prop == "pressure":
              ROUT = os.path.join(outpath,"Vapor_Pressure.list")
           self.setenv("ROUT",ROUT)
           self.setenv("AVERAGE_TYPE","LJ")
           self.setenv("AVERAGE_TYPE","LJ")
           res = os.system("R --slave < %s" %(gro_R))
           if res != 0:
              self.errorexit()
           pressure_file = os.path.join(tmppath,"average.txt")
           pressure_vector = self.read_last_parameter(pressure_file)
           pressure = self.mean(pressure_vector)
           try:
                if prop == "P":
		   os.stat(os.path.join(outpath,"Pressure.sd"))
                else:
                   os.stat(os.path.join(outpath,"Vapor_Pressure.sd"))
           except:
                dp = self.sd(pressure_vector)
           else:
                if prop == "P":
		   sd_file = os.path.join(outpath,"Pressure.sd")
                else:
                   sd_file = os.path.join(outpath,"Vapor_Pressure.sd")
		dp = self.read_last_parameter(sd_file)[0]
       
           property = pressure
           noise = dp

        # Diffusion Coefficient (general or for cation (IL))
        if prop == "diff_cat":
           DC_file_cat = os.path.join(outpath,"D_cat.erg")
           DC_list_cat = self.read_line(DC_file_cat,listing=True)
           DC_cat = float(DC_list_cat[0])
           try:
		os.stat(os.path.join(outpath,"D_cat.sd"))
           except:
                noise = 0.1*DC_cat
           else:
                sd_file = os.path.join(outpath,"D_cat.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = DC_cat

        # Diffusion Coefficient for anion (IL)
        if prop == "diff_an":          
           DC_file_an = os.path.join(outpath,"D_an.erg")
           DC_list_an = self.read_line(DC_file_an,listing=True)
           DC_an = float(DC_list_an[0])
           try:
		os.stat(os.path.join(outpath,"D_an.sd"))
           except:
                noise = 0.1*DC_an
           else:
                sd_file = os.path.join(outpath,"D_an.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = DC_an

        # Reorientation Time
        if prop == "reorient_time":
           REOR_file = os.path.join(outpath,"tau2.erg")
           REOR_list = self.read_line(REOR_file,listing=True)
           REOR = float(REOR_list[0])
           try:
		os.stat(os.path.join(outpath,"tau2.sd"))
           except:
                noise = 0.1*REOR
           else:
                sd_file = os.path.join(outpath,"tau2.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = REOR
        
        # Viscosity
        if prop == "viscosity":
           VISC_file = os.path.join(outpath,"visc.erg")
           VISC_list = self.read_line(VISC_file,listing=True)
           VISC = float(VISC_list[0])
           try:
		os.stat(os.path.join(outpath,"visc.sd"))
           except:
                noise = 0.1*VISC
           else:
                sd_file = os.path.join(outpath,"visc.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = VISC

        # Thermal conductivity
        if prop == "thermal_cond":
           Tcond_file = os.path.join(outpath,"tcond.erg")
           Tcond_list = self.read_line(Tcond_file,listing=True)
           Tcond = float(Tcond_list[0])
           try:
		os.stat(os.path.join(outpath,"tcond.sd"))
           except:
                noise = 0.1*Tcond
           else:
                sd_file = os.path.join(outpath,"tcond.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = Tcond

        # Isothermal compressibility
        if prop == "isothermal_compress":
           Icomp_file = os.path.join(outpath,"icomp.erg")
           Icomp_list = self.read_line(Icomp_file,listing=True)
           Icomp = float(Icomp_list[0])
           try:
		os.stat(os.path.join(outpath,"icomp.sd"))
           except:
                noise = 0.1*Icomp
           else:
                sd_file = os.path.join(outpath,"icomp.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = Icomp

        # Volume expansivity
        if prop == "volume_expans":
           Vexp_file = os.path.join(outpath,"vexp.erg")
           Vexp_list = self.read_line(Vexp_file,listing=True)
           Vexp = float(Vexp_list[0])
           try:
		os.stat(os.path.join(outpath,"vexp.sd"))
           except:
                noise = 0.1*Vexp
           else:
                sd_file = os.path.join(outpath,"vexp.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = Vexp

        # dH/dP
        if prop == "dHdP":
           DHdP_file = os.path.join(outpath,"dHdP.erg")
           DHdP_list = self.read_line(DHdP_file,listing=True)
           DHdP = float(DHdP_list[0])
           try:
		os.stat(os.path.join(outpath,"dHdP.sd"))
           except:
                noise = 0.1*DHdP
           else:
                sd_file = os.path.join(outpath,"dHdP.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = DHdP

        # speed of sound
        if prop == "speed":
           Speed_file = os.path.join(outpath,"speed.erg")
           Speed_list = self.read_line(Speed_file,listing=True)
           Speed = float(Speed_list[0])
           try:
		os.stat(os.path.join(outpath,"speed.sd"))
           except:
                noise = 0.1*Speed
           else:
                sd_file = os.path.join(outpath,"speed.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = Speed

        # isobaric heat capacity
        if prop == "isobaric_heatcap":
           IbarHeat_file = os.path.join(outpath,"ibarHeat.erg")
           IbarHeat_list = self.read_line(IbarHeat_file,listing=True)
           IbarHeat = float(IbarHeat_list[0])
           try:
		os.stat(os.path.join(outpath,"ibarHeat.sd"))
           except:
                noise = 0.1*IbarHeat
           else:
                sd_file = os.path.join(outpath,"ibarHeat.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = IbarHeat
        
        # isochoric heat capacity
        if prop == "isochoric_heatcap":
           IchorHeat_file = os.path.join(outpath,"ichorHeat.erg")
           IchorHeat_list = self.read_line(IchorHeat_file,listing=True)
           IchorHeat = float(IchorHeat_list[0])
           try:
		os.stat(os.path.join(outpath,"ichorHeat.sd"))
           except:
                noise = 0.1*IchorHeat
           else:
                sd_file = os.path.join(outpath,"ichorHeat.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = IchorHeat

        # dU/dV
        if prop == "dUdV":
           DUdV_file = os.path.join(outpath,"dUdV.erg")
           DUdV_list = self.read_line(DUdV_file,listing=True)
           DUdV = float(DUdV_list[0])
           try:
		os.stat(os.path.join(outpath,"dUdV.sd"))
           except:
                noise = 0.1*DUdV
           else:
                sd_file = os.path.join(outpath,"dUdV.sd")
		noise = self.read_last_parameter(sd_file)[0]

           property = DUdV

        return [property,noise]

    def get_properties_from_files(self,properties,name,P=None):
        """ collocates all required properties from MD tool output files """

        properties_list = []
        noise_list = []

        program = os.getenv("program")

        T_total = []
        if os.getenv("VLE") == "y":
           VLE_T_range = eval(os.getenv("VLETRANGE"))
           for k in range(len(VLE_T_range)):
               T_total.append(VLE_T_range[k])
        if os.getenv("transport") == "TRUE":
           T_trans = eval(os.getenv("TTRANS"))
           for k in range(len(T_trans)):
               T_total.append(T_trans[k])
        if os.getenv("BULKNPT") == "y":
           T_bulk_NPT = eval(os.getenv("NPTTRANGE"))
           for k in range(len(T_bulk_NPT)):
               T_total.append(T_bulk_NPT[k])
        if os.getenv("BULKNVT") == "y":
           T_bulk_NVT = eval(os.getenv("NVTTRANGE"))
           for k in range(len(T_bulk_NVT)):
               T_total.append(T_bulk_NVT[k])

        property_dict = eval(os.getenv("PROPERTYDICT"))
        for k in range(len(T_total)):
            if P == None:
               outpath = os.path.join(os.getenv("TMPPATH"),repr(T_total[k]),name)
            else:
               outpath = os.path.join(os.getenv("TMPPATH"),repr(P),repr(T_total[k]),name)
            if program == "towhee" and "vapor" in properties:
               NMOL_file = os.path.join(outpath,"liquid_molecules.txt")
               N = int(self.read_line(NMOL_file))
            else:
               N = None
            properties_isothermal = []
            noise_isothermal = []
            prop_list = property_dict[repr(T_total[k])]
            for prop in prop_list:
                [property,noise] = self.make_properties_file_by_lists(outpath,prop)
                properties_isothermal.append(property)
                noise_isothermal.append(noise)
            properties_list.append(properties_isothermal)
            noise_list.append(noise_isothermal)

        return [properties_list,noise_list]

    def get_density_from_files(self,P_range,T_range,properties,properties_list,noise_list,name,P):
        """ collocates all densities at all pressures and temperatures from MD tool output files """

        rho_list = []
        drho_list = []
        for p in range(len(P_range)):
            rho_list_isobaric = []
            drho_list_isobaric = []
            if P_range[p] != P:
               for k in range(len(T_range)):
                   outpath = os.path.join(os.getenv("TMPPATH"),repr(P_range[p]),repr(T_range[k]),name)
                   [rho,drho] = self.make_properties_file_by_lists(outpath,"density")
                   rho_list_isobaric.append(rho)
                   drho_list_isobaric.append(drho)
            else:
                   rho_index = properties.index("density")
                   for k in range(len(T_range)):
                       rho_list_isobaric.append(properties_list[k][rho_index])
                       drho_list_isobaric.append(noise_list[k][rho_index])
            rho_list.append(rho_list_isobaric)
            drho_list.append(drho_list_isobaric)

        return [rho_list,drho_list]

    def get_properties_file(self,program,name=None):
        """ produces the file containing the simulated properties and returns its name """

        eff = os.getenv("efficiency")

        if program == "korr":
           properties_file = os.path.join(os.getenv("OUTPATH"),"properties_cp.txt")
        else:
            properties_file = os.path.join(os.getenv("OUTPATH"),"properties.txt")
        if program != "korr":  
            print "Creating properties files..."       
            if eff == "y":
                 if name.startswith("original") or name.startswith("armijo"):
                    self.setenv("fits",None)
                 else:
                    self.setenv("fits","y")
            if eff == None or (eff == "y" and (name.startswith("original") or name.startswith("armijo"))):
              properties = eval(os.getenv("PROPERTIES"))
              if os.getenv("density_fit") == None or os.getenv("density_fit") == "linear":
                [properties_list,noise_list] = self.get_properties_from_files(properties,name)
                rho_list = []
                drho_list = []
                
              if os.getenv("fits") == "y":                               
                 eff_obj = efficiency.Efficiency()
                 eff_obj.fits(properties_list,noise_list,rho_list,drho_list,name,self)         
    
              outpath = os.getenv("ORIG_OUTPATH")
              properties_file = os.path.join(outpath,"properties.txt")
            else: # effiziente Gradienten-/Hesseberechnung mit red. Parametern
              T_range = eval(os.getenv("TRANGE_ERROR"))
              properties = eval(os.getenv("PROPERTIES"))
              eff_obj = efficiency.Efficiency()
              i = len(self.epsilon_ref_tilde_list)-1
              outpath = os.getenv("ORIG_OUTPATH")
              red_properties_file = os.path.join(outpath,"reduced_properties.txt")
              if os.getenv("density_fit") == None or os.getenv("density_fit") == "linear":
                 [properties_list,noise_list] = eff_obj.compute_properties(red_properties_file,self.epsilon_ref_tilde_list[i],self.sigma_ref_tilde_list[i],T_range,properties)
                 rho_list = []
                 drho_list = []
                 P_range_actual = None
                 
              T_range_actual = eff_obj.actual_T_range(T_range,self.epsilon_ref,self.epsilon_ref_tilde_list[i])

              eff_obj.fits(properties_list,noise_list,rho_list,drho_list,name,self,T_range_actual,P_range_actual)

              outpath = os.getenv("ORIG_OUTPATH")
              properties_file = os.path.join(outpath,"properties.txt")
                
            if os.getenv("fits") == None:
               property_dict = eval(os.getenv("PROPERTYDICT"))
               T_total = []
               if os.getenv("VLE") == "y":
                  VLE_T_range = eval(os.getenv("VLETRANGE"))
                  for k in range(len(VLE_T_range)):
                      T_total.append(VLE_T_range[k])
               if os.getenv("transport") == "TRUE":
                  T_trans = eval(os.getenv("TTRANS"))
                  for k in range(len(T_trans)):
                      T_total.append(T_trans[k])
               if os.getenv("BULK_NPT") == "y":
                  T_bulk_NPT = eval(os.getenv("NPTTRANGE"))
                  for k in range(len(T_bulk_NPT)):
                      T_total.append(T_bulk_NPT[k])
               if os.getenv("BULK_NVT") == "y":
                  T_bulk_NVT = eval(os.getenv("NVTTRANGE"))
                  for k in range(len(T_bulk_NVT)):
                      T_total.append(T_bulk_NVT[k])
               f = open(properties_file,"w")
               for k in range(len(T_total)):
                   for l in range(len(property_dict[repr(T_total[k])])):
                       f.write(repr(properties_list[k][l])+"\n")
               f.close()              

        if os.getenv("METHOD") not in ["stoll","stoll2"]:
           return properties_file
        else:  
           outpath = os.getenv("ORIG_OUTPATH")
           properties_range_file = os.path.join(outpath,"properties_range.txt")
           return [properties_file,properties_range_file]

    def append_parameter(self,x,filename,rep=True):
        """ appends parameter x to file """

        f = open(filename,"a")

        for i in range(len(x)):
            if i != len(x)-1:
               if rep == True:
                  f.write(repr(x[i])+" ")
               else:
                  f.write(x[i]+" ")
            else:
               if rep == True:
                  f.write(repr(x[i])+"\n")
               else:
                  f.write(x[i]+"\n")

        f.close()
    
        return 0

    def floatable(self,prop):
        """ makes properties floatable (sometimes , is used instead of .) """
        
        fprop = ""

        for i in range(len(prop)):
            if prop[i] != ",":
               fprop += prop[i]
            else:
               fprop += "."

        return fprop

    def get_properties(self, properties_file):
        """ reads physical properties from a file """

        act_prop = []
        
        f = open(properties_file,"r")
        prop = f.readlines()
        f.close()

        

        for i in range(len(prop)):
            prop[i] = prop[i][:(len(prop[i])-1)]
            prop[i] = self.floatable(prop[i])
            act_prop.append(float(prop[i]))

    
        return act_prop

    def loss_function(self,properties_file):
        """ calculates the loss function """

        properties = self.get_properties(properties_file)
        targets = self.get_properties(os.getenv("TARGET"))
        weights = eval(os.getenv("WEIGHTS"))

        loss = 0.0

        for i in range(len(properties)):
            calc = properties[i]
            exp = targets[i]
            #if calc/exp < 0.1 or calc/exp > 10:
             #  print calc
             #  print exp
             #  print "Property too far away from experiment. Please check your previous simulation run!"
              # loss = 100000
            summand = pow(1.0-(float(calc)/float(exp)),2)
            loss += weights[i]*summand

        #loss = math.sqrt(loss)

        eps = float(os.getenv("LIMIT"))

        if loss <= eps:
           print "Loss function less than %f." %(eps)

        return loss

    def loss_function_2(self,properties_file,weight_file):
        """ calculates a different loss function """

        properties = self.get_properties(properties_file)
        targets = self.get_properties(os.getenv("TARGET"))
        weights = self.get_properties(weight_file)

        loss = 0.0

        for i in range(len(properties)):
            calc = properties[i]
            exp = targets[i]
            w = weights[i]
            summand = (1.0/(w*exp))*pow(calc - exp,2)
            loss += summand

        loss = loss/len(properties)

        if loss <= eps:
           print "Loss function less than %f." %(eps)

        return loss

    def evaluation(self,gradient,Hessian):
        """ evaluates the optimization """

        epsilon = float(os.getenv("LIMIT"))

        norm_gradient = self.norm(gradient)
 
        if norm_gradient <= epsilon:
           return "ok"
        else:
           return "no"

        # todo: Hessian positve definit ???
           
    def which_parameter(self,order,parind_list,k):
        """ finds out if the actual parameter is epsilon, sigma, or a charge """


        for i in range(len(parind_list)):
            if k in parind_list[i]:
               par = order[i]

        return par

    def get_all_kind_of(self,par_name,x,order,parind_list):
        """ gets all parameters of kind par_name in vector x """

        par_list = []
        for i in range(len(x)):
            par_list.append(self.which_parameter(order,parind_list,i+1))

        kind_list = []
        index_list = []
        for i in range(len(par_list)):
            if par_list[i] == par_name:
               kind_list.append(x[i])
               index_list.append(i)
 
        return [kind_list,index_list]

    def index_in_list_of_lists(self,elem,ll):
        """ finds out where an element is situated within a list of lists """

        for i in range(len(ll)):
            if elem in ll[i]:
               return i

        return False

    def change_parameter(self,y,h,i):
        """ changes parameter x to x+h at index i and writes new parameter into parameter file """

        order = eval(os.getenv("PARAMORDER"))
        number = eval(os.getenv("PARAMNUMBER"))
        parind_list = eval(os.getenv("PARINDLIST"))


        if os.getenv("efficiency") == None:
           y[i] = y[i] + h
           self.append_parameter(y,os.getenv("PARAMFILE"))
           if os.getenv("CHARGEOPT") == "True":
              rest_charge = self.get_rest_charges(y)
              rest_charge_file = os.getenv("PARAMFILE")+"_rest_charges"          
              self.append_parameter(rest_charge,rest_charge_file)

        else:
           eff_obj = efficiency.Efficiency()
           if i == 0:
              self.epsilon_ref_tilde_list = []
              self.sigma_ref_tilde_list = []
           res = eff_obj.change_parameter(y,h,i,self.par_red_list)
           if res[0] == "epsilon" or res[0] == "charge":
              self.epsilon_ref_tilde_list.append(res[1])
              self.sigma_ref_tilde_list.append(self.sigma_ref)
           if res[0] == "sigma":
              self.sigma_ref_tilde_list.append(res[1])
              self.epsilon_ref_tilde_list.append(self.epsilon_ref)


        if os.getenv("METHOD") == "stoll" or os.getenv("RETURN") == "y":
           return y

        else:
           return

    def change_two_parameters(self,y,h,i,j):
        """ changes x to x+h at indices i and j and writes new parameter into parameter file """
  
        order = eval(os.getenv("PARAMORDER"))
        number = eval(os.getenv("PARAMNUMBER"))
        parind_list = eval(os.getenv("PARINDLIST"))

        if os.getenv("efficiency") == None:
           y[i] = y[i] + h
           y[j] = y[j] + h
           self.append_parameter(y,os.getenv("PARAMFILE"))
           if os.getenv("CHARGEOPT") == "True":
              rest_charge = self.get_rest_charges(y)
              rest_charge_file = os.getenv("PARAMFILE")+"_rest_charges"          
              self.append_parameter(rest_charge,rest_charge_file)
        else:
           eff_obj = efficiency.Efficiency()
           res = eff_obj.change_two_parameters(y,h,i,j,self.par_red_list)
           if (res[0] == "epsilon" and res[1] == "epsilon") or (res[0] == "charge" and res[1] == "charge"):
              self.epsilon_ref_tilde_list.append(res[2])
              self.sigma_ref_tilde_list.append(self.sigma_ref)
           if res[0] == "sigma" and res[1] == "sigma":
              self.sigma_ref_tilde_list.append(res[2])
              self.epsilon_ref_tilde_list.append(self.epsilon_ref)
           if res[0] == "sigma" and (res[1] == "epsilon" or res[1] == "charge"):
              self.sigma_ref_tilde_list.append(res[2])
              self.epsilon_ref_tilde_list.append(res[3])
           if (res[0] == "epsilon" or res[0] == "charge") and res[1] == "sigma":
              self.sigma_ref_tilde_list.append(res[3])
              self.epsilon_ref_tilde_list.append(res[2])
           if (res[0] == "epsilon" and res[1] == "charge") or (res[0] == "charge" and res[1] == "epsilon"):
              self.epsilon_ref_tilde_list.append(res[2])
              self.sigma_ref_tilde_list.append(self.sigma_ref)
 
        return
       
    def merge_properties(self, properties_file_list,d):
        """ merges several properties files to one (d = number of properties) """ 


        properties_matrix = []
        for i in range(d):
            properties_matrix.append([])
            for j in range(len(properties_file_list)):
                properties_matrix[i].append(0.0)
            
        for i in range(len(properties_file_list)):
            properties_file = properties_file_list[i]
            properties = self.get_properties(properties_file)
            for j in range(d):
                prop = properties[j]
                properties_matrix[j][i] = prop

        tmppath = os.getenv("TMPPATH")
        properties_var_file = os.path.join(tmppath,"properties_var.txt")
        
        self.write_matrix(properties_matrix,properties_var_file)

        return properties_var_file


    def read_matrix(self,matrix_file,fl=True):
        """ reads a matrix from file """

        f = open(matrix_file,"r")

        matrix = []

        while True:
              line = f.readline()
              if not line:
                 break
              line = string.split(line," ")
              matrix.append(line)
        f.close()

        for i in range(len(matrix)):
            x = matrix[i]
            x[len(x)-1] = x[len(x)-1][:(len(x[len(x)-1])-1)]
            if fl == True:
               matrix[i] = self.float_vector(x)

        return matrix

    def write_vector(self,vector,filename):
        """ writes vector into file """

        f = open(filename,"w")

        for i in range(len(vector)):
            f.write(repr(vector[i])+"\n")

        f.close()

    def write_matrix(self,matrix,filename):
        """ writes matrix into file """

        f = open(filename,"w")
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if j != len(matrix[i])-1:
                   f.write(repr(matrix[i][j])+" ")
                else:
		   f.write(repr(matrix[i][j])+"\n")
        f.close()

        return 0

    def quadratic_expression(self,v,A):
        """ calculates the quadratic expression (vT)Av of a vector v and a matrix A """
 
        qe = 0.0
 
        if len(v) != len(A):
           print "Number of rows must coincide with the dimension of the vector:"
           print "Matrix A has %d rows but vector v is of dimension %d." %(len(A),len(v))
           self.errorexit()

        if len(A[0]) != len(v):
           print "Number of columns must coincide with the dimension of the vector:"
           print "Matrix A has %d columns but vector v is of dimension %d." %(len(A[0]),len(v))
           self.errorexit()

        for i in range(len(A)):
            for j in range(len(v)):
                qe += A[i][j]*v[i]*v[j]

        return qe

    def scalar_product(self,v1,v2):
        """ calculates the scalar product of two vector v1, v2 """
 
        if len(v1) != len(v2):
           print "Error in %s scalar_product: Vectors do not have the same length." %("opt_class.py")
           self.errorexit()
 
        prod = 0.0

        for i in range(len(v1)):
            prod += v1[i]*v2[i]

        return prod

    def mean(self,v):
        """ calculates the average of a float vector v
        """

        m = 0.0
        for i in range(len(v)):
            m += v[i]
    
        return m/len(v)

    def sd(self,v):
        """ calculates the standard deviation of a sample v 
        """

        sd = 0.0
   
        m = self.mean(v)
     
        for i in range(len(v)):
            sd += pow(v[i]-m,2)
 
        if len(v) > 1:
           sd = sd/(len(v)-1)

        return math.sqrt(sd)

    def norm(self,v):
        """ calculates euclidean norm of vector v """

        norm_v = 0
    
        for i in range(len(v)):
            norm_v += pow(v[i],2)

        norm_v = math.sqrt(norm_v)
    
        return(norm_v)

    def L2_norm(self,int_err):
        """ computes the L2 norm of an interpolation error """

        L2 = 0.0

        for i in range(len(int_err)):
            L2 += pow(int_err[i],2)

        return math.sqrt(L2)

    def abs_vector(self,vector):
        """ returns [|v1|,|v2|,...,|vN|], where vi are the components of a vector """

        vector2 = []

        for i in range(len(vector)):
            vector2.append(abs(vector[i]))
     
        return vector2

    def Linf_norm(self,int_err):
        """ computes the Linf norm of an interpolation error """

        return max(self.abs_vector(int_err))

    def frobenius_norm(self,A):
        """ computes the Frobenius norm of a matrix A """

        fn = 0.0

        for i in range(len(A)):
            for j in range(len(A[i])):
                fn += pow(A[i][j],2)

        return math.sqrt(fn)

    def faculty(self,x):
        """ faculty of x """

        if x > 1:
           return x * self.faculty(x-1)
        else:
           return 1

    def binom(self,n,k):
        """ binomial coefficient """

        return self.faculty(n)/(self.faculty(k)*self.faculty(n-k))

    def vector_sum(self,v,w):
        """ computes the sum of two vectors v and w """

        if len(v) != len(w):
           print "Error while computing the sum of 2 vectors: v and w do not have the same length."
           print "v has length %d and w has length %d." %(len(v),len(w))
           self.errorexit()

        s = []

        for i in range(len(v)):
            s.append(v[i]+w[i])

        return s

    def matrix_sum(self,A,B):
        """ computes the sum of two matrices A and B """

        if len(A) != len(B):
           print "Error while computing the sum of 2 matrices: A and B do not have the same number of rows."
           print "A: %d and B: %d." %(len(A),len(B))
           self.errorexit()

        if len(A[0]) != len(B[0]):
           print "Error while computing the sum of 2 matrices: A and B do not have the same number of columns."
           print "A: %d and B: %d." %(len(A[0]),len(B[0]))
           self.errorexit()

        C = []

        for i in range(len(A)):
            C.append([])
            for j in range(len(A[i])):
                C[i].append(A[i][j]+B[i][j])

        return C

    def mult_vector(self,l,v):
        """ multiplies a vector v by l """

        lv = []
 
        for i in range(len(v)):
            lv.append(l*v[i])

        return lv

    def mult_matrix(self,l,A):
        """ multiplies a matrix A by l """

        for i in range(len(A)):
            for j in range(len(A[i])):
                A[i][j] = l*A[i][j]

        return A

    def matrix_vector(self,A,v):
        """ computes the matrix vector product of matrix A and vector v """

        r = []

        n = len(A)
        m = len(A[0])

        for i in range(n):
            s = 0.0
            for j in range(m):
                s += A[i][j]*v[j]
            r.append(s)

        return r


    def complete_summary_table(self,loop,properties,targets,loss):
        """ completes the summary table by actual information """

        table_file = os.path.join(os.getenv("GROW_HOME"),"table.txt")

        if loop == 0:
           try:
               os.stat(table_file)
           except:
               pass
           else:
               os.system("rm -f %s" %(table_file))

        print "Completing table file %s" %(table_file)

        if loop == 0:
           sum_f = open(table_file,"w")
        else:
           sum_f = open(table_file,"a")

        row = "%d " %(loop)
         
        for i in range(len(properties)):
            if i != len(properties)-1:
               row += "%.4f %.4f %.4f %.4f " %(properties[i],targets[i],properties[i]-targets[i],(properties[i]-targets[i])/targets[i])
            else:
               row += "%.4f %.4f %.4f %.4f %.4f\n" %(properties[i],targets[i],properties[i]-targets[i],(properties[i]-targets[i])/targets[i],loss)
        sum_f.write(row)

        sum_f.close()

        return

    def complete_summary_file(self,loop,param,gradient,norm,properties,targets,loss):
        """ completes the summary file by actual information """

        method = os.getenv("METHOD")

        summary_file = os.path.join(os.getenv("GROW_HOME"),"summary.txt")

        if loop == 0:
           try:
               os.stat(summary_file)
           except:
               pass
           else:
               os.system("rm -f %s" %(summary_file))

        print "Completing summary file %s" %(summary_file)

        # get temperatures from working directories
        T_range_err = []
        fl = False
        working_directories = self.read_matrix(os.path.join(os.getenv("GROW_HOME"),"working_directories.txt"), fl)
        for wd_list in working_directories:
            working_directory = wd_list[0]
            temperature = float(os.path.basename(working_directory))
            if temperature not in T_range_err:
               T_range_err.append(temperature)

        sum_f = open(summary_file,"a")
        sum_f.write("\n")
        sum_f.write("Values for x%d:\n" %(loop))
        sum_f.write("x%d = " %(loop))
        for i in range(len(param)):
            if i != len(param)-1:
               sum_f.write(repr(param[i])+" ")
            else:
               sum_f.write(repr(param[i])+"\n")

        if method not in ["dgafo","simplex","EA"]:
           sum_f.write("gradient = ")
           for i in range(len(gradient)):
               if i != len(gradient)-1:
                  sum_f.write(repr(gradient[i])+" ")
               else:
                  sum_f.write(repr(gradient[i])+"\n")
           sum_f.write("norm of gradient = %f\n" %(norm))

        if method == "dgafo":
           sum_f.write("Interpolated function value: "+repr(gradient)+"\n")
           sum_f.write("Delta: "+repr(norm)+"\n")  
           if loss != 0.0:
              sum_f.write("Relative SG Interpolation Error: %f\n" %((gradient-loss)/loss))      

        # get property dict
        property_dict = eval(self.read_line(os.path.join(os.getenv("GROW_HOME"),"property_dict.txt")))
        MAPE_dict = {}
        if properties != []:

           i = 0
           for T in T_range_err:

               prop_names = property_dict[repr(T)]

               for prop_name in prop_names:
                   if prop_name not in MAPE_dict.keys():
                      MAPE_dict[prop_name] = []

               for prop_name in prop_names:
                   sum_f.write("T = %.1f: calc %s = %.4f, exp %s = %.4f, abs error = %.4f, rel error = %.4f\n" %(T,prop_name,properties[i],prop_name,targets[i],properties[i]-targets[i],(properties[i]-targets[i])/targets[i]))
                   prop_index = prop_names.index(prop_name)
                   MAPE_dict[prop_name].append(abs((targets[i]-properties[i])/targets[i])*100)
                   i += 1
        
           for prop_name in MAPE_dict:
               MAPE_list = list(MAPE_dict[prop_name])
               sum_f.write("MAPE on %s: %.4f\n" %(prop_name,sum(MAPE_list)/float(len(MAPE_list))))

        sum_f.write("loss = %f\n" %(loss))
        sum_f.close()

        return

    def write_summary_file(self,x,f_x,loop,nof_param):  
        """ writes the results into the summary file """

        print "See details in %s." %(os.path.join(os.getenv("OPTOUTPATH"),"summary.txt"))
        sum_f = open(os.path.join(os.getenv("OPTOUTPATH"),"summary.txt"),"a")
        sum_f.write("\n")
        sum_f.write("Optimal set of parameters:")
        for i in range(len(x)):
            sum_f.write(" %.4f" %(x[i]))
        sum_f.write("\n")
        sum_f.write("Value of loss function: %.4f\n" %(f_x))
        sum_f.write("Number of iterations: %d\n" %(loop))
        sum_f.write("Number of function evaluations: %d\n" %(nof_param))
        sum_f.close()

        return 0

    def test_independence(self,base,loop):
        """ tests base vectors for linear independence """

        tmppath = os.getenv("TMPPATH")
        A_file = os.path.join(tmppath,"Base.txt")
        self.write_matrix(base,A_file)

        independence_script = os.getenv("INDEPSCRIPT")

        res = os.system("R --slave < %s" %(independence_script))
        if res != 0:
           self.errorexit()
        
        evaluation = self.read_line(os.path.join(tmppath,"R_eval.txt"))

        if evaluation == "dep":
           print "The base vectors for the efficient gradient calculation are not linearly independent. Please modify the optimization parameter 'h' by one percent or so and set the restart parameter to '%d 1'." %(loop)
           self.errorexit()

        return

    def matrix_column(self,j,matrix):
        """ return column j of a matrix """

        column = []

        for i in range(len(matrix)):
            column.append(matrix[i][j])

        return column

    def transpose(self,A):
        """ transposes the matrix A """

        AT = []

        for i in range(len(A)):
            AT.append([])
            for j in range(len(A[i])):
                AT[i].append(A[j][i])

        return AT

    def base_transformation(self,vector,o,base,solvescript=None):
        """ computes a base transformation of a cartesian vector into a base with origin o """

        tmppath = os.getenv("TMPPATH")
        if solvescript == None:
           solvescript = os.getenv("SOLVESCRIPT")

        A = self.transpose(base)

        b = self.vector_sum(vector,self.mult_vector(-1.0,o))

        self.write_matrix(A,os.path.join(tmppath,"A.txt")) 
        self.write_vector(b,os.path.join(tmppath,"b.txt"))
        
        res = os.system("R --slave < %s" %(solvescript))
        if res != 0:
           self.errorexit()

        trans_vector = self.read_last_parameter(os.path.join(tmppath,"v.txt"))

        return trans_vector

    def cartesian_coordinates(self,vector,o,base):
        """ performs a transformation of a vector w.r.t. a base with origin o into cartesian coordinates """

        cart_vector = []

        for i in range(len(o)):
            S = 0.0
            for j in range(len(base)):
                S += vector[j]*base[j][i]
            S += o[i]
            cart_vector.append(S)

        return cart_vector

    def evaluate_quadratic_form(self,v,c,g,H):
        """ calculates c + gTv + 1/2vTHv """

        return c + self.scalar_product(g,v) + 0.5*self.quadratic_expression(v,H)


    def get_KKT_matrix(self,H,l):
        """ computes KKT matrix H + 2lI """

        KKT_matrix = []
        for i in range(len(H)):
            KKT_matrix.append(list(H[i]))

        for i in range(len(KKT_matrix)):
            KKT_matrix[i][i] = KKT_matrix[i][i]+2*l

        return KKT_matrix
        
    def reduce_quadratic_form(self,d,Delta,g,H,R_pd_script):
        """ reduces the quadratic form subject to norm(d) <= Delta following a TR theorem """
  

        sp = self.scalar_product(g,d)
        norm_d = self.norm(d)

        if sp > 0.0:
           new_d = self.mult_vector(-Delta/norm_d,d)
 

        if sp <= 0.0:
           self.setenv("ev_method","neg_ev")
           res = os.system("R --slave < %s" %(R_pd_script))
           if res != 0:
              self.errorexit()
           candidate_z = self.read_last_parameter(os.path.join(os.getenv("TMPPATH"),"eigenvector.txt"))
           if self.scalar_product(g,candidate_z) <= 0:
              z = list(candidate_z)
           else:
              z = self.mult_vector(-1,candidate_z)
           norm_z = self.norm(z)
           
           if norm_d < Delta:
              p = 2*self.scalar_product(z,d)/pow(norm_z,2)
              q = (pow(norm_d,2)-pow(Delta,2))/pow(norm_z,2)
              alpha = max(-p/2 + math.sqrt( pow(p/2,2) - q ),-p/2 - math.sqrt( pow(p/2,2) - q ))
              new_d = self.vector_sum(d,self.mult_vector(alpha,z))
           else:
              if self.scalar_product(z,d) != 0.0 and False:
                 new_d = self.vector_sum(d,self.mult_vector(-2*self.scalar_product(z,d)/pow(norm_z,2),z))
              else:
                 qe = self.quadratic_expression(z,H)
                 p = 2*self.scalar_product(g,z)/abs(qe)
                 q = -abs(self.scalar_product(g,z))/abs(qe)
                 alpha = -p/2 + math.sqrt( (pow(p/2,2) - q) ) + 1.0
                 dalpha = self.vector_sum(d,self.mult_vector(alpha,z))
                 new_d = self.vector_sum(d,self.mult_vector(-2*pow(Delta,2)/(pow(Delta,2)+pow(alpha,2)*pow(norm_z,2)),dalpha))
          
        return new_d

    def comp_dev_sim(self,loop,restart,dim):
        """ finds out from the restart parameter, which component deviations should be simulated """

        if restart == [0]:
           sim_indic = 1
        else:
           if loop < restart[0]:
              sim_indic = dim+1
           elif loop == restart[0]:
              if len(restart) == 2:
                 sim_indic = restart[1]
              elif len(restart) > 2:
                 sim_indic = dim+1
              else:
                 sim_indic = 1
           else: 
              sim_indic = 1

        return sim_indic

    def tr_sim(self,tr,loop,restart,dim):
        """ finds out from the restart parameter, when the solution of the TR subproblem should be simulated """
        
        if restart == [0]:
           sim_indic = True
        else:
           if loop < restart[0]:
              sim_indic = False
           elif loop == restart[0]:
              if len(restart) == 5:
                 if restart[4] == 0:
                    sim_indic = True
                 else:
                    if tr < restart[4]:
                       sim_indic = False
                    else:
                       sim_indic = True
              elif len(restart) == 3:
                 if restart[2] == 0:
                    sim_indic = True
                 else:
                    if tr < restart[2]:
                       sim_indic = False
                    else:
                       sim_indic = True
              else:
                 sim_indic = True
           else:
              sim_indic = True

        return sim_indic

    def comp_dev_Hess_diag_sim(self,loop,restart,dim):
        """ finds out from the restart parameter, which component deviations for the diagonal entries of the Hessian should be simulated """

        if restart == [0]:
           sim_indic = 1
        else:
           if loop < restart[0]:
              sim_indic = dim+1
           elif loop == restart[0]:
              if len(restart) == 4:
                 if restart[2] == restart[3]:
                    sim_indic = restart[2]
                 else:
                    sim_indic = dim+1
              elif len(restart) >= 5:
                 sim_indic = dim+1
              else:
                 sim_indic = 1
           else: 
              sim_indic = 1

        return sim_indic

    def comp_dev_Hess_nondiag_sim(self,loop,restart,dim):
        """ finds out from the restart parameter, which component deviations for the nondiagonal entries of the Hessian should be simulated """

        if restart == [0]:
           sim_indic = [1,1]
        else:
           if loop < restart[0]:
              sim_indic = [dim+1,dim+1]
           elif loop == restart[0]:
              if len(restart) == 4:
                 if restart[2] != restart[3]:
                    sim_indic = [restart[2],restart[3]]
                 else:
                    sim_indic = [1,1]
              elif len(restart) >= 5:
                 sim_indic = [dim+1,dim+1]
              else:
                 sim_indic = [1,1]
           else: 
              sim_indic = [1,1]

        return sim_indic

    def gradient_calc(self,loop,restart):
        """ finds out from the restart parameter, if the gradient should be calculated or not """

        if restart == [0]:
           gradient_indic = True
        else:
           if loop < restart[0]:
              gradient_indic = False
           elif loop == restart[0]:
              if len(restart) == 2:
                 gradient_indic = True
              elif len(restart) > 2:
                 if restart[2] == 0:
                    gradient_indic = True
                 else:
                    gradient_indic = False
              else:
                 gradient_indic = True
           else:
              gradient_indic = True

        return gradient_indic

    def Hessian_calc(self,loop,restart):
        """ finds out from the restart parameter, if the Hessian should be calculated or not """

        if restart == [0]:
           Hessian_indic = True
        else:
           if loop < restart[0]:
              Hessian_indic = False
           elif loop == restart[0]:
              if len(restart) == 4:
                 Hessian_indic = True
              elif len(restart) > 4:
                 if restart[4] == 0:
                    Hessian_indic = True
                 else:
                    Hessian_indic = False
              else:
                 Hessian_indic = True
           else:
              Hessian_indic = True

        return Hessian_indic

    def delete(self):
       """ deletes temporary files, if required """

       print "Deleting temporary files, if required."
       d = os.getenv("delete")

       if d == "dt":
          os.removedirs(os.getenv("TMPPATH"))
       if d == "dte":
          os.removedirs(os.getenv("OUTPATH"))
       if d == "d":
          os.removedirs(os.getenv("TMPPATH"))
          os.removedirs(os.getenv("OUTPATH"))
 
          return

    def __del__(self):
          """
          destructor
          """	
	  del self	
