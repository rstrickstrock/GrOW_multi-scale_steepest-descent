#     class Trace, 2014/08/21 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		trace.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		21-08-2014
#
#     Description:	error handling functionality
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
#			traceback
#			time
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
import sys
import traceback
import time

# ***********************************************************************
# class Trace
# ***********************************************************************

class Trace:
       
    def __init__(self):
          """ constructor """
          
	  self.__trace_obj = None


    def errorexit(self,tool=False):
       """gives out general error message in case of severe error"""

       if tool == False:
          print "Severe error, aborting."
       else:
          print "Simulation tool error, aborting" 

       if tool == False:
          del self

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

    def __del__(self):
          """
          destructor
          """	
	  del self
