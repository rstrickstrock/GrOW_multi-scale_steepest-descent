#     class System, 2014/08/21 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		system.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		21-08-2014
#
#     Description:	system and environment control functionality
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

# ***********************************************************************
# class System
# ***********************************************************************

class System:
       
    def __init__(self):
          """ constructor """
          
	  self.__system_obj = None

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

    def __del__(self):
          """
          destructor
          """	
	  del self

