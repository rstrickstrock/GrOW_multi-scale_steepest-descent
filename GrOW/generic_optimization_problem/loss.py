#     class Loss, 2014/04/30 huelsmann FhI SCAI/HBRS
#
#     start_doc
#
#     Script:		loss.py
#  
#     Author:		Marco Huelsmann
# 
#     Date:		30-04-2014
#
#     Description:	class Loss (derived from Objective_Function)
#
#     Usage:		by defining an instance
#
#     Arguments:	
#
#     Options:	       
#			
#     Output:	        
#
#     Imported:		
#
#     Calling:	
#			
#
#     Modifications:
#
#	$Log: md_opt$
#
#     end_doc

# ***********************************************************************
# required python modules
# ***********************************************************************

from objective_function import Objective_Function

# ***********************************************************************
# class Loss
# ***********************************************************************

class Loss(Objective_Function):
       
    def __init__(self, dimension, weights, targets):
          """ constructor """
          
          self._dimension = dimension
          self._gradients_of_estimations = []
          self._weights = weights
          self._targets = targets

          if self._objective_function == "Physical_Properties_Loss" or "QM_MM_Loss":
             self._weights = weights
             self._targets = targets
             #print "self._weights:",self._weights
             #print "self._targets:",self._targets
             self._estimations = []

             self._hessians_of_estimations = []

          if self._objective_function == "PhysProp_QMMM_Loss":
             #self._weights_QM_MM = [weights[0]]
             #print "loss.py: self._weights_QM_MM:",self._weights_QM_MM
             #self._weights_Physical_Properties = [weights[1]]
             #print "loss.py: self._weights_Physical_Properties:",self._weights_Physical_Properties
             #print "loss.py: len(self._weights_Physical_Properties):",len(str(self._weights_Physical_Properties))
             #self._weights = self._weights_QM_MM
             #for w in self._weights_Physical_Properties:
             #   self._weights.append(w)
            # print "loss.py: self._weights:",self._weights

             #self._targets_QM_MM = targets[0]
             #print "loss.py: self._targets_QM_MM:",self._targets_QM_MM
             #self._targets_Physical_Properties = targets[1]
             #print "loss.py: self._targets_Physical_Properties:",self._targets_Physical_Properties
             #self._targets = self._targets_QM_MM
             #for t in self._targets_Physical_Properties:
             #   self._targets.append(t)
            # print "loss.py: self._targets:",self._targets

             self._estimations_QM_MM = []
             self._estimations_Physical_Properties = []
             ################### check if distinguishing these vars is necessary ###################
             self._hessians_of_estimations_QM_MM = []
             self._hessians_of_estimations_Physical_Properties = []

          #tmp_file=open('/home/rstric2s/targets.txt','w')
          #tmp_file.write("targets")
          #for i in range(len(self._targets)):
          #   tmp_file.write(str(self._targets[i]) + "\n")
          #tmp_file.close

          #tmp_file=open('/home/rstric2s/weights.txt','w')
          #tmp_file.write("weights")
          #for i in range(len(self._weights)):
          #   tmp_file.write(str(self._weights[i]) + "\n")
          #tmp_file.close

          #Objective_Function.__init__(self, dimension)


    def get_dimension(self):
        """ return the dimension of the loss function """
        return self._dimension

    def check_target(self, target, i, epsilon=0.000001):
        """ checks if target is below epsilon (close to zero) to prevent dividing by zero. """
        #target = 0.00000000001
        #print "set target to %s" %(repr(target))
        if target < epsilon:
            print "In loss.check_target: To prevent a division by zero or numerical instability target[%s] = %s has been changed to %s" %(repr(i),repr(target),repr(epsilon))
            target = epsilon

        return target

    def get_weights(self):
        """ return the weights of the loss function """
        return self._weights

    def get_weights_QM_MM(self):
        """ return QM - weights of the loss function """
        #print "loss.py: return - self._weights_QM_MM:",self._weights_QM_MM
        return self._weights_QM_MM

    def get_weights_Physical_Properties(self):
        """ return Physical Properties - weights of the loss function """
        #print "loss.py: return - self.weights_Physical_Properties:",self._weights_Physical_Properties
        return self._weights_Physical_Properties

    def get_estimations(self):
        """ return the estimations of the loss function """
        return self._estimations

    def get_estimation_gradients(self):
        """ returns gradients of estimations """
        return self._gradients_of_estimations

    def get_estimation_gradient(self, i):
        """ returns gradient of estimation i """
        return self._gradients_of_estimations[i]

    def get_estimation_derivative(self, i, k):
        """ returns kth derivative of estimation i """
        return self._gradients_of_estimations[i][k]

    def get_estimation_hessians(self):
        """ returns hessians of estimations """
        return self._hessians_of_estimations

    def get_estimation_hessian(self, i):
        """ returns hessian of estimation i """
        return self._hessians_of_estimations[i]

    def get_estimation_hessian_row(self, i, k):
        """ returns kth row of hessian of estimation i """
        return self._hessians_of_estimations[i][k]

    def get_estimation_hessian_entry(self, i, k, l):
        """ returns entry kl of hessian of estimation i """
        return self._hessians_of_estimations[i][k][l]

    def set_estimations(self, x):
        """ sets estimations, overloaded method """
        pass

    def set_estimation_gradients(self, x):
        """ sets gradients of estimations, overloaded method """
        pass

    def set_estimation_hessians(self, x):
        """ sets hessians of estimations, overloaded method """
        pass

    def get_targets(self):
        """ return the targets of the loss function """
        return self._targets
    
    def get_function_values(self, parameter_set):
        """ returns a list of function values for a given parameter set"""
        pass
    
    def get_function_value(self, x):
        """ return the loss function value """

        loss_value = 0.0

        #print "in loss.py: self._estimations:",self._estimations
        #print "in loss.py: self._targets:",self._targets

        for i in range(len(self._targets)):
            #print "estimations[%i]: %s" %(i, self._estimations[i]) 
            #print "targets[%i]: %s" %(i, self._targets[i]) 
            summand = pow(1.0-(float(self._estimations[i])/float(self._targets[i])),2)
            loss_value += self._weights[i]*summand

        return loss_value

    def get_gradient(self, x):
        """ return the gradient of the loss function """
         
        gradient = []
        #print "in loss.py: self._targets =",self._targets
        #print "in loss.py: len(self._targets) =",str(len(self._targets))
        #print "in loss.py: self._weights =",self._weights
        #print "in loss.py: len(self._weights) =",str(len(self._weights))
        #test_file = open('/home/rstric2s/properties_loss.py.txt','w')
        
        for k in range(self._dimension):
            #test_file.write("k \t i \t weight[i] \t goe \t target[i] \t estimation[i]\n")
            partial_F = 0.0
            for i in range(len(self._targets)):
                #tmp_prop = str(k) + "\t" + str(i) + "\t" + str(self._weights[i]) + "\t" + str(self._gradients_of_estimations[i][k]) + "\t" + str(self._targets[i]) + "\t" +  str(self._estimations[i]) + "\n"
                #test_file.write(tmp_prop)
                self._targets[i] = self.check_target(self._targets[i], i)
                #print "in loss.py: self._weights[i]:", self._weights[i], type(self._weights[i])
                #print "in loss.py: self._gradients_of_estimations[i][k]:", self._gradients_of_estimations[i][k], type(self._gradients_of_estimations[i][k])
                #print "in loss.py: self._targets[i]:", self._targets[i], type(self._targets[i])
                #print "in loss.py: self._estimations[i]:", self._estimations[i], type(self._estimations[i])
                partial_F += self._weights[i] * (self._gradients_of_estimations[i][k] * ( (self._targets[i]-self._estimations[i])/pow(self._targets[i],2) ) )
            partial_F = - 2*partial_F
            gradient.append(partial_F)
            #test_file.write("\n")
        #test_file.close()

        return gradient

    def get_hessian(self, x):
        """ return the Hessian of the loss function """
        
        hessian = []

        # initialization
        for k in range(self._dimension):
            hessian.append([])
            for l in range(self._dimension):
                hessian[k].append(l)

        # calculation
        for k in range(self._dimension):
            for l in range(self._dimension):
                for i in range(len(self._targets)):
                    hessian[k][l] += self._weights[i]/pow(self._targets[i],2) * ( 4 * self._gradients_of_estimations[i][l] - 2 * self._hessians_of_estimations[i][k][l] * (self._targets[i] - self._estimations[i]) )

        return hessian

    def get_boundary(self):
	""" get boundary vectors """
	pass

    def __del__(self):
          """
          destructor
          """	
	  del self	
