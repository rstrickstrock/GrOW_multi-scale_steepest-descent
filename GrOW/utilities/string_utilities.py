#     class String_Utilities, 2014/08/21 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		string_utilities.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		21-08-2014
#
#     Description:	string handling functionality
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

import os
import os.path
import string
import sys

# ***********************************************************************
# class String_Utilities
# ***********************************************************************

class String_Utilities:
       
    def __init__(self):
          """ constructor """
          
	  self.__string_utilities_obj = None

    def eval_vector(self,vector,check=False):
        """ converts a vector to a string vector, if possible """
 
        for i in range(len(vector)):
            if check == False:
               vector[i] = eval(vector[i])
            else:
               	try:
   		   vector[i] = eval(vector[i])
		except:
		   return vector[i]

        return vector

    def floatable(self,string_value):
        """ makes a string value floatable (sometimes , is used instead of .) """
        
        floatable_value = ""

        for i in range(len(string_value)):
            if string_value[i] != ",":
               floatable_value += string_value[i]
            else:
               floatable_value += "."

        return floatable_value

    def float_vector(self,vector,check=False):
        """ converts a vector to a float vector, if possible """
 
        for i in range(len(vector)):
            if check == False:
               vector[i] = float(vector[i])
            else:
               	try:
   			vector[i] = float(vector[i])
		except:
			return vector[i]
        return vector

    def __del__(self):
          """
          destructor
          """	
	  del self
