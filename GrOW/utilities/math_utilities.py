#     class Math, 2014/08/21 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		math.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		21-08-2014
#
#     Description:	mathematical functionality
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
# 			
#     Output:	        
#
#     Imported:		math
#			os
#			random
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
import random as rndm
import sys

from utilities.io import IO
i_o = IO()

from utilities.system import System
sy = System()

from utilities.trace import Trace
tr = Trace()

# ***********************************************************************
# class Math
# ***********************************************************************

class Math:
       
    def __init__(self):
          """ constructor """
          
	  self.__math_obj = None

    def quadratic_expression(self,v,A):
        """ calculates the quadratic expression (vT)Av of a vector v and a matrix A """
 
        qe = 0.0
 
        if len(v) != len(A):
           print "Number of rows must coincide with the dimension of the vector:"
           print "Matrix A has %d rows but vector v is of dimension %d." %(len(A),len(v))
           tr.errorexit()

        if len(A[0]) != len(v):
           print "Number of columns must coincide with the dimension of the vector:"
           print "Matrix A has %d columns but vector v is of dimension %d." %(len(A[0]),len(v))
           tr.errorexit()

        for i in range(len(A)):
            for j in range(len(v)):
                qe += A[i][j]*v[i]*v[j]

        return qe

    def scalar_product(self,v1,v2):
        """ calculates the scalar product of two vector v1, v2 """
 
        if len(v1) != len(v2):
           print "Error in %s scalar_product: Vectors do not have the same length." %("opt_class.py")
           tr.errorexit()
 
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
           tr.errorexit()

        s = []

        for i in range(len(v)):
            s.append(v[i]+w[i])

        return s

    def matrix_sum(self,A,B):
        """ computes the sum of two matrices A and B """

        if len(A) != len(B):
           print "Error while computing the sum of 2 matrices: A and B do not have the same number of rows."
           print "A: %d and B: %d." %(len(A),len(B))
           tr.errorexit()

        if len(A[0]) != len(B[0]):
           print "Error while computing the sum of 2 matrices: A and B do not have the same number of columns."
           print "A: %d and B: %d." %(len(A[0]),len(B[0]))
           tr.errorexit()

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

    def test_independence(self,base,loop):
        """ tests base vectors for linear independence """

        A_file = os.path.join(os.getenv("GROW_HOME"),"Base.txt")
        i_o.write_matrix(base,A_file)

        independence_script = os.getenv("INDEPSCRIPT")

        res = os.system("R --slave < %s" %(independence_script))
        if res != 0:
           tr.errorexit()
        
        evaluation = self.read_line(os.path.join(os.getenv("GROW_HOME"),"R_eval.txt"))

        if evaluation == "dep":
           print "The base vectors for the efficient gradient calculation are not linearly independent. Please modify the optimization parameter 'h' by one percent or so and set the restart parameter to '%d 1'." %(loop)
           tr.errorexit()

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

        if solvescript == None:
           solvescript = os.getenv("SOLVESCRIPT")

        A = self.transpose(base)

        b = self.vector_sum(vector,self.mult_vector(-1.0,o))

        i_o.write_matrix(A,os.path.join(os.getenv("GROW_HOME"),"A.txt")) 
        i_o.write_vector(b,os.path.join(os.getenv("GROW_HOME"),"b.txt"))
        
        res = os.system("R --slave < %s" %(solvescript))
        if res != 0:
           tr.errorexit()

        trans_vector = self.read_last_parameter(os.path.join(os.getenv("GROW_HOME"),"v.txt"))

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
           sy.setenv("ev_method","neg_ev")
           res = os.system("R --slave < %s" %(R_pd_script))
           if res != 0:
              tr.errorexit()
           candidate_z = i_o.read_last_parameter(os.path.join(os.getenv("GROW_HOME"),"eigenvector.txt"))
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

    def index_in_list_of_lists(self,elem,ll):
        """ finds out where an element is situated within a list of lists """

        for i in range(len(ll)):
            if elem in ll[i]:
               return i

        return False

    def introduce_uncertainties(self,value,uncertainty):
        """ introduces artificial uncertainty on a value """

        dev_percent = rndm.uniform(0.0,uncertainty)

        dev = rndm.choice([-1,1])*dev_percent*value
        
        value_biased = value + dev

        return value_biased

    def __del__(self):
          """
          destructor
          """	
	  del self
