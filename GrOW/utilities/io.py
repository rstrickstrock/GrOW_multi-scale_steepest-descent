#     class IO, 2014/08/21 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		io.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		21-08-2014
#
#     Description:	file I/O functionality
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
# 			
#     Output:	        
#
#     Imported:		os
#			os.path
#			string
#			sys
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

import math
import os
import string
import sys

from utilities.string_utilities import String_Utilities
su = String_Utilities()

from utilities.trace import Trace
tr = Trace()

# ***********************************************************************
# class IO
# ***********************************************************************

class IO:
       
    def __init__(self):
          """ constructor """
          
	  self.__io_obj = None

    def read_config_file(self,config_file):
       """
       reads configuration file "config_file" and returns a config object
       """

       try:
	  import ConfigParser
       except:
	  print "No module named ConfigParser but required."
	  print "Abort."
	  tr.errorexit()

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
	      tr.errorexit()

       return source_file_list

    def read_dim(self):
       """ reads the dimension from the parameter file """

       param_file = os.getenv("PARAMFILE")

       x = os.popen("tail -1 %s" %(param_file)).readlines()
       x = string.split(x[0]," ")

       return len(x)

    def read_last_parameter(self,filename,fl=True):
        """ reads the last parameter from a file """


        x = os.popen("tail -1 %s" %(filename)).readlines()
        x = string.split(x[0]," ")
        x[len(x)-1] = x[len(x)-1][:(len(x[len(x)-1])-1)]
        if fl == True:
           x = su.float_vector(x)

        return x

    def read_last_two_parameters(self,filename,fl=True):
        """ reads the last two parameters from a file """

        dim = int(os.getenv("dim"))

        x = os.popen("tail -2 %s" %(filename)).readlines()
        x1 = string.split(x[0]," ")
        x2 = string.split(x[1]," ")
        x1[len(x1)-1] = x1[len(x1)-1][:(len(x1[len(x1)-1])-1)]
        if fl == True:
           x1 = su.float_vector(x1)
        x2[len(x2)-1] = x2[len(x2)-1][:(len(x2[len(x2)-1])-1)]
        if fl == True:
           x2 = su.float_vector(x2)

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
    
        return

    def get_properties(self, properties_file):
        """ reads physical properties from a file """

        act_prop = []
        
        f = open(properties_file,"r")
        prop = f.readlines()
        f.close()

        for i in range(len(prop)):
            prop[i] = prop[i][:(len(prop[i])-1)]
            prop[i] = su.floatable(prop[i])
            act_prop.append(float(prop[i]))

    
        return act_prop

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
            # eliminate blanks
            y = []
            for k in range(len(x)):
                if x[k] != "":
                   y.append(x[k])
            if fl == True:
               matrix[i] = su.float_vector(y)

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

    def complete_summary_table(self,loop,properties,targets,loss):
        """ completes the summary table by actual information """

        table_file = os.path.join(os.getenv("GROW_HOME"),"table.txt")

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

    def complete_summary_file(self,loop,param,gradient,norm,properties,targets,loss,hybrid_opt=False):
        """ completes the summary file by actual information """

        method = os.getenv("METHOD")

        sum_f = open(os.path.join(os.getenv("GROW_HOME"),"summary.txt"),"a")
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

        if properties != []:
           T_range_err = os.getenv("TRANGE_ERROR")
           if T_range_err == None:
              T_range_err = [float(os.getenv("TEMPERATURE"))]
           else:
              T_range_err = eval(T_range_err)

           prop_names = eval(os.getenv("PROPERTIES"))

           MAPE_list = []

           for prop_name in prop_names:
               MAPE = []
               for T in T_range_err:
                   MAPE.append(0.0)
               MAPE_list.append(MAPE)

           i=0
           for T in T_range_err:
               T_index = T_range_err.index(T)
               for prop_name in prop_names:
                   sum_f.write("T = %.1f: calc %s = %.4f, exp %s = %.4f, abs error = %.4f, rel error = %.4f\n" %(T,prop_name,properties[i],prop_name,targets[i],properties[i]-targets[i],(properties[i]-targets[i])/targets[i]))
                   prop_index = prop_names.index(prop_name)
                   MAPE_list[prop_index][T_index] = abs((targets[i]-properties[i])/targets[i])*100
                   i += 1
        
           for prop_name in prop_names:
               prop_index = prop_names.index(prop_name)
               sum_f.write("MAPE on %s: %.4f\n" %(prop_name,self.mean(MAPE_list[prop_index])))

        sum_f.write("loss = %f\n" %(loss))
        sum_f.close()

        return

    def write_summary_file(self,x,f_x,loop,nof_param):  
        """ writes the results into the summary file """

        print "See details in %s." %(os.path.join(os.getenv("OPTOUTPATH"),"summary.txt"))
        sum_f = open(os.path.join(os.getenv("GROW_HOME"),"summary.txt"),"a")
        sum_f.write("\n")
        sum_f.write("Optimal set of parameters:")
        for i in range(len(x)):
            sum_f.write(" %.4f" %(x[i]))
        sum_f.write("\n")
        sum_f.write("Value of loss function: %.4f\n" %(f_x))
        sum_f.write("Number of iterations: %d\n" %(loop))
        sum_f.write("Number of function evaluations: %d\n" %(nof_param))
        sum_f.close()

        return

    def __del__(self):
          """
          destructor
          """	
	  del self
