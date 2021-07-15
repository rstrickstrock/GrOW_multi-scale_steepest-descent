import os
import os.path
import sys
import traceback
import string
import math
import opt_class


class Efficiency:
       
    def __init__(self):
          """ constructor """
          self.obj = None
          self.opt = opt_class.Optimization()


    def rho0(self,T,R):
        """ help function for density fit """

        rho_0 = R[0] + R[1]*T + R[2]*pow(T,2) + R[3]*pow(T,3)

        return rho_0

    def BT(self,T,B):
        """ help function for density fit """

        B_T = B[0] + B[1]*T + B[2]*pow(T,2)

        return B_T

    def calculate_density(self,T_range,rho_coefs,P=None):
        """ calculates the density for different temperatures from the density nl regression coefficients """

        rho_list = []

        P0 = 1.0
        if P == None:
           P = float(os.getenv("PRESSURE"))
        A = rho_coefs[0]
        B = [rho_coefs[1],rho_coefs[2],rho_coefs[3]]
        R = [rho_coefs[4],rho_coefs[5],rho_coefs[6],rho_coefs[7]]

        for T in T_range:
            rho0 = self.rho0(T,R)
            B_T = self.BT(T,B)
            rho_list.append(rho0/(1-A*math.log((B_T+P)/(B_T+P0))))


        return rho_list

    def calculate_linear_density(self,T_range_err,rho_coefs):
        """ calculates linearly fitted density """

        rho_list = []

        for i in range(len(T_range_err)):
            rho_list.append(rho_coefs[0]+rho_coefs[1]*T_range_err[i])

        return rho_list

    def evaluate_density(self,rho_coefs,T_range,P_range,rho_list,w):
        """ calculates the density error function for different temperatures and pressures from the actual density nl regression coefficients """

        err = 0.0       

        calc_list = []
        for P in P_range:
            calc_list.append(self.calculate_density(T_range,rho_coefs,P))

        for p in range(len(P_range)):
            for k in range(len(T_range)):
                #err += pow((rho_list[p][k]-calc_list[p][k])/(w[p][k]*rho_list[p][k]),2)
                err += pow((rho_list[p][k]-calc_list[p][k])/w[p][k],2)

        return err

    def calc_deriv_rho_A(self,rho_coefs,T_range,P_range):
        """ partial rho/partial A """

        deriv_rho_A = []

        P0 = 1.0
        A = rho_coefs[0]
        B = [rho_coefs[1],rho_coefs[2],rho_coefs[3]]
        R = [rho_coefs[4],rho_coefs[5],rho_coefs[6],rho_coefs[7]]

        for P in P_range:
            deriv_rho_A_isobaric = []
            for T in T_range:
                rho_0 = self.rho0(T,R)
                B_T = self.BT(T,B)
                external = rho_0/pow((1-A*math.log((B_T+P)/(B_T+P0))),2)
                internal = math.log((B_T+P)/(B_T+P0))
                deriv_rho_A_isobaric.append(external*internal)
            deriv_rho_A.append(deriv_rho_A_isobaric)

        return deriv_rho_A

    def calc_deriv_rho_B(self,rho_coefs,T_range,P_range,j):
        """ partial rho/partial Bj """

        deriv_rho_B = []

        P0 = 1.0
        A = rho_coefs[0]
        B = [rho_coefs[1],rho_coefs[2],rho_coefs[3]]
        R = [rho_coefs[4],rho_coefs[5],rho_coefs[6],rho_coefs[7]]

        for P in P_range:
            deriv_rho_B_isobaric = []
            for T in T_range:
                rho_0 = self.rho0(T,R)
                B_T = self.BT(T,B)
                external = rho_0/pow((1-A*math.log((B_T+P)/(B_T+P0))),2)
                internal = (A*pow(T,j)*(P0-P))/((B_T+P)*(B_T+P0))
                deriv_rho_B_isobaric.append(external*internal)
            deriv_rho_B.append(deriv_rho_B_isobaric)

        return deriv_rho_B

    def calc_deriv_rho_rho0(self,rho_coefs,T_range,P_range,j):
        """ partial rho/partial rho0j """

        deriv_rho_rho0 = []

        P0 = 1.0
        A = rho_coefs[0]
        B = [rho_coefs[1],rho_coefs[2],rho_coefs[3]]

        for P in P_range:
            deriv_rho_rho0_isobaric = []
            for T in T_range:
                B_T = self.BT(T,B)
                deriv = pow(T,j)/(1-A*math.log((B_T+P)/(B_T+P0)))
                deriv_rho_rho0_isobaric.append(deriv)
            deriv_rho_rho0.append(deriv_rho_rho0_isobaric)

        return deriv_rho_rho0


    def calculate_error_gradient(self,x,T_range,P_range,rho_list,w):
        """ calculates the error gradient dependent on the actual density nl regression coefficients """

        grad = []
        
        deriv_rho_A = self.calc_deriv_rho_A(x,T_range,P_range)

        deriv_rho_B = []
        for j in range(3):
            deriv_rho_B.append(self.calc_deriv_rho_B(x,T_range,P_range,j))

        deriv_rho_rho0 = []
        for j in range(4):
            deriv_rho_rho0.append(self.calc_deriv_rho_rho0(x,T_range,P_range,j))

        grad_rho = [deriv_rho_A]

        for j in range(3):
            grad_rho.append(deriv_rho_B[j])

        for j in range(4):
            grad_rho.append(deriv_rho_rho0[j])

        calc_list = []
        for P in P_range:
            calc_list.append(self.calculate_density(T_range,x,P))

        for i in range(8):
            grad_comp = 0.0
            for p in range(len(P_range)):
                for k in range(len(T_range)):
                    grad_comp += w[p][k]*(rho_list[p][k]-calc_list[p][k])*grad_rho[i][p][k]
            grad.append(-2.0*grad_comp)

        return grad

    def armijo(self,x_old,err,grad,d,beta,sigma,T_range,P_range,rho_list,w):
        """ performs Armijo step length control for SD optimization for density fit """

        l = 1       
        while True:
              sl = pow(beta,l)
              x = self.opt.vector_sum(x_old,self.opt.mult_vector(sl,d))

              err_sl = self.evaluate_density(x,T_range,P_range,rho_list,w)

              if err_sl <= err+sigma*sl*self.opt.scalar_product(grad,d):
                 break
              l += 1



        return x

    def density_fit(self,T_range,P_range,rho_list,drho_list):
        """ calculates the density nl regression coefficients """
       
        print "Density fit..."

        w = []
       
        for i in range(len(rho_list)):
            wi = []
            for j in range(len(rho_list[i])):
                wi.append(1.0/(drho_list[i][j]))
                #wi.append(drho_list[i][j])
                #wi.append(1.0)
            w.append(wi)


        P0 = 0.1
	#x = [8.2858*pow(10,-2),7.6867*pow(10,2),-2.7198,3.2094*pow(10,-3),1.6909*pow(10,-3),-2.9885,7.0444*pow(10,-3),-7.1575*pow(10,-6)] B in MPa, rho in kg/m3
        x = [8.2858*pow(10,-2),7.6867*pow(10,3),-20.7198,3.2094*pow(10,-2),1.6909*pow(10,3),-2.9885,7.0444*pow(10,-3),-7.1575*pow(10,-6)]
      
        #for i in range(7): 	# 7 = Anzahl Koeff.
        #    x.append(1.0)

        eps = 0.001
        err_old = -1.0

        print "Solving nonlinear optimization problem by SD method ... "

        k = 0
        while True: 
              err = self.evaluate_density(x,T_range,P_range,rho_list,w)
              grad = self.calculate_error_gradient(x,T_range,P_range,rho_list,w)
              if self.opt.norm(grad) <= eps or abs(err-err_old) <= eps:
                 break
              else:
                 err_old = err
              d = self.opt.mult_vector(-1.0,grad)
              if self.opt.norm(grad) >= 1.0:
                 d = self.opt.mult_vector(1.0/self.opt.norm(d),d)
              beta = 0.01
              sigma = 0.0
              x = self.armijo(x,err,grad,d,beta,sigma,T_range,P_range,rho_list,w)
              k += 1            

        return x

    def linear_density_fit(self,T_range,P_range,rho_list,drho_list,name):
        """ performs a linear fit for the density """

        print "Linear Density fit..."

        tmppath = os.getenv("TMPPATH")
        path = os.path.join(tmppath,"density_fits",name)
        
        try:
		os.stat(path)
        except:
		os.makedirs(path)

        w = []

        for i in range(len(drho_list)):
            w.append(1.0/drho_list[i])
            #w.append(1.0)

        A = []

        for T in T_range:
            A.append(T)

        b = []
 
        for rho in rho_list:
            b.append(rho)

        f1 = open(os.path.join(path,"A.txt"),"w")
        for i in range(len(A)):
            f1.write(repr(A[i])+"\n")
        f1.close()
            
        f2 = open(os.path.join(path,"w.txt"),"w")
        for i in range(len(w)):
            f2.write(repr(w[i])+"\n")
        f2.close()    
 
        f3 = open(os.path.join(path,"b.txt"),"w")
        for i in range(len(b)):
            f3.write(repr(b[i])+"\n")
        f3.close()  

        del f1
        del f2
        del f3

        # fit by weighted linear regression
        fit_script = os.getenv("FITSCRIPT")

        os.system("R --slave < %s --args %s" %(fit_script,path))

        rho_coefs = self.opt.read_last_parameter(os.path.join(path,"coefs.txt"))

        os.system("rm %s %s %s" %(os.path.join(path,"A.txt"),os.path.join(path,"w.txt"),os.path.join(path,"b.txt")))

        return rho_coefs

        

    def guggenheim_fit(self,T_range,sld_list,dsld_list,path,epsilon_ref,sigma_ref):
        """ performs a Guggenheim fit for the reduced saturated liquid density """

        print "Guggenheim fit..."

        for i in range(len(sld_list)):
            sld_list[i] = self.reduce_density(sld_list[i],sigma_ref)

        for i in range(len(T_range)):
            T_range[i] = self.reduce_temperature(T_range[i],epsilon_ref)

        path = os.path.join(path,"guggenheim_fits")
        
        try:
		os.stat(path)
        except:
		os.makedirs(path)

        w = []
       
        for i in range(len(sld_list)):
            w.append(1.0/dsld_list[i])
            #w.append(1.0)

        coefs_list = []
        residuals_list = []

        for p in range(49):
            Tc = (1+(p+1)*0.01)*max(T_range)

            Tr_list = []
            for i in range(len(T_range)):
                Tr_list.append(T_range[i]/Tc)

            A = []
            b = []

            # A
            for i in range(len(T_range)):
                T = T_range[i]
                row = [pow(Tc-T,(1.0/3.0)),Tc-T,pow(Tc-T,1.5)]
                A.append(row)
 
            # b
            b = list(sld_list)

            self.opt.write_matrix(A,os.path.join(path,"A.txt"))
            f1 = open(os.path.join(path,"b.txt"),"w")
            for i in range(len(b)):
                f1.write(repr(b[i])+"\n")
            f1.close()
            f2 = open(os.path.join(path,"w.txt"),"w")
            for i in range(len(w)):
                f2.write(repr(w[i])+"\n")
            f2.close()    
            del f1
            del f2

            # guggenheim fit by weighted linear regression
            fit_script = os.getenv("FITSCRIPT")

            os.system("R --slave < %s --args %s" %(fit_script,path))

            sld_coefs = self.opt.read_last_parameter(os.path.join(path,"coefs.txt"))
            sld_residuals = self.opt.read_last_parameter(os.path.join(path,"residuals.txt"))
            coefs_list.append(sld_coefs)
            residuals_list.append(self.opt.norm(sld_residuals))

        min_ind = residuals_list.index(min(residuals_list))

        Tc = (1+(min_ind+1)*0.01)*max(T_range)
        Tc_phys = self.temperature(Tc,epsilon_ref)
        sld_coefs = coefs_list[min_ind]

        fTc = open(os.path.join(path,"critical_temperatures.txt"),"w")
        fTc.write(repr(Tc)+repr(Tc_phys)+"\n")
        fTc.close()

        f3 = open(os.path.join(path,"sld_coefs.txt"),"w")
        for i in range(len(sld_coefs)):
            if i != len(sld_coefs)-1:
               f3.write(repr(sld_coefs[i])+" ")
            else:
               f3.write(repr(sld_coefs[i])+"\n")
        f3.close()
        del f3

        os.system("rm %s %s %s" %(os.path.join(path,"A.txt"),os.path.join(path,"w.txt"),os.path.join(path,"b.txt")))

        return [sld_coefs,Tc]

    def pressure_fit(self,T_range,p_list,dp_list,path):
        """ calculates the vapor pressure nl regression coefficients """

        print "Vapor pressure fit..."        
        
        try:
		os.stat(path)
        except:
		os.makedirs(path)

        w = []
       
        for i in range(len(p_list)):
            w.append(1.0/dp_list[i])
            #w.append(1.0)        

        A = []

        for T in T_range:
            row = [1.0/T,1.0/pow(T,4)]
            A.append(row)

        b = []

        for p in p_list:
            b.append(math.log(p))

        self.opt.write_matrix(A,os.path.join(path,"A.txt"))
            
        f2 = open(os.path.join(path,"w.txt"),"w")
        for i in range(len(w)):
            f2.write(repr(w[i])+"\n")
        f2.close()    
 
        f3 = open(os.path.join(path,"b.txt"),"w")
        for i in range(len(b)):
            f3.write(repr(b[i])+"\n")
        f3.close()  

        del f2
        del f3

        # fit by weighted linear regression
        fit_script = os.getenv("FITSCRIPT")

        os.system("R --slave < %s --args %s" %(fit_script,path))

        p_coefs = self.opt.read_last_parameter(os.path.join(path,"coefs.txt"))

        os.system("rm %s %s %s" %(os.path.join(path,"A.txt"),os.path.join(path,"w.txt"),os.path.join(path,"b.txt")))

        return p_coefs
    

    def vapor_fit(self,T_range,vapor_list,dvapor_list,path):
        """ calculates the enthalpy of vaporization nl regression coefficients """
       
        print "Enthalpy of vaporization fit..."
        
        try:
		os.stat(path)
        except:
		os.makedirs(path)

        w = []
       
        for i in range(len(vapor_list)):
            w.append(1.0/dvapor_list[i])
            #w.append(1.0)            

        coefs_list = []
        residuals_list = []

        for p in range(49):
            Tc = (1+(p+1)*0.01)*max(T_range)

            Tr_list = []
            for i in range(len(T_range)):
                Tr_list.append(T_range[i]/Tc)

            A = []
            b = []

            # A
            for i in range(len(Tr_list)):
                Tr = Tr_list[i]
                #row = [math.log(1-Tr),-Tr]
                row = [math.log(1-Tr),math.log(1-Tr)*Tr,math.log(1-Tr)*pow(Tr,2)]
                A.append(row)
 
            # b
            for i in range(len(vapor_list)):
                b.append(math.log(vapor_list[i]))

            self.opt.write_matrix(A,os.path.join(path,"A.txt"))
            f1 = open(os.path.join(path,"b.txt"),"w")
            for i in range(len(b)):
                f1.write(repr(b[i])+"\n")
            f1.close()
            f2 = open(os.path.join(path,"w.txt"),"w")
            for i in range(len(w)):
                f2.write(repr(w[i])+"\n")
            f2.close()    
            del f1
            del f2

            # fit by weighted linear regression
            fit_script = os.getenv("FITSCRIPT")

            os.system("R --slave < %s --args %s" %(fit_script,path))
 
            vapor_coefs = self.opt.read_last_parameter(os.path.join(path,"coefs.txt"))
            vapor_residuals = self.opt.read_last_parameter(os.path.join(path,"residuals.txt"))
            coefs_list.append(vapor_coefs)
            residuals_list.append(self.opt.norm(vapor_residuals))

        min_ind = residuals_list.index(min(residuals_list))

        Tc = (1+(min_ind+1)*0.01)*max(T_range)
        vapor_coefs = coefs_list[min_ind]

        fTc = open(os.path.join(path,"critical_temperature.txt"),"w")
        fTc.write(repr(Tc)+"\n")
        fTc.close()

        os.system("rm %s %s %s" %(os.path.join(path,"A.txt"),os.path.join(path,"w.txt"),os.path.join(path,"b.txt")))

        return [vapor_coefs,Tc]

    def calculate_sld(self,T_range,Tc,sld_coefs,epsilon_ref,sigma_ref):
        """ computes the saturated liquid density from Guggenheim coefficients """

        sld_list = []

        for i in range(len(T_range)):
            T = self.reduce_temperature(T_range[i],epsilon_ref)
            sld_list.append(sld_coefs[0]+sld_coefs[1]*pow(Tc-T,(1.0/3.0))+sld_coefs[2]*(Tc-T)+sld_coefs[3]*pow(Tc-T,1.5))

        for i in range(len(sld_list)):
            sld_list[i] = self.density(sld_list[i],sigma_ref)

        return sld_list

    def calculate_reduced_sld(self,T_range,Tc,sld_coefs,epsilon_ref):
        """ computes the saturated liquid density from Guggenheim coefficients """

        sld_list = []
        T_range_red = []

        for i in range(len(T_range)):
            T = self.reduce_temperature(T_range[i],epsilon_ref)
            T_range_red.append(T)
            sld_list.append(sld_coefs[0]+sld_coefs[1]*pow(Tc-T,(1.0/3.0))+sld_coefs[2]*(Tc-T)+sld_coefs[3]*pow(Tc-T,1.5))

        return [sld_list,T_range_red]

    def calculate_pressure(self,T_range,p_coefs):
        """ computes the fitted vapor pressure """

        p_list = []

        alpha = p_coefs[0]
        beta = p_coefs[1]
        gamma = p_coefs[2]

        for T in T_range:
            p_list.append(math.exp(alpha+beta/T+gamma/pow(T,4)))

        return p_list

    def calculate_vapor(self,T_range,Tc,vapor_coefs):
        """ computes the fitted enthalpy of vaporization """
      
        vapor_list = []

        A = math.exp(vapor_coefs[0])
        beta = vapor_coefs[1]
        alpha = vapor_coefs[2]
        gamma = vapor_coefs[3] #

        for i in range(len(T_range)):
            Tr = T_range[i]/Tc
            #vapor_list.append(A*pow(1-Tr,beta)*math.exp(-alpha*Tr))
            vapor_list.append(A*pow(1-Tr,beta+alpha*Tr+gamma*pow(Tr,2)))

        return vapor_list


    def calc_T_range_for_plotting(self,T_range_actual):
        """ calculates the temperature range for the fit plots """

        tmppath = os.getenv("TMPPATH")
	T_file = os.path.join(tmppath,"plot_temperature.txt")

        try:
		os.stat(T_file)
        except:
                lower_T_bound = min(T_range_actual)
                upper_T_bound = max(T_range_actual)

       		T_range = []
    
                nof_T = int(upper_T_bound-lower_T_bound+1)

                for i in range(nof_T):
                    T_range.append(lower_T_bound+i)

                fT = open(T_file,"w")
                for i in range(len(T_range)):
                    if i != len(T_range)-1:
                       fT.write(repr(T_range[i])+" ")
                    else:
                       fT.write(repr(T_range[i])+"\n")
                fT.close()
        else:
		T_range = self.opt.read_last_parameter(T_file)

        return T_range

    def plot_fits(self,coefs,T_range,T_range_actual,property_list,noise_list,name,old_property_list,property,Tc=None,epsilon_ref=None,sigma_ref=None):
        """ makes a plot of the desired fit """

        print "Plotting the fits..."

        tmppath = os.getenv("TMPPATH")
        plot_fits_script = os.getenv("PLOTFITSSCRIPTS")

        if property == "density":
           if os.getenv("density_fit") == "specific":
              property_range = self.calculate_density(T_range,coefs)
           else:
              property_range = self.calculate_linear_density(T_range,coefs)
           path = os.path.join(tmppath,"density_fits",name)
        if property == "vapor":
           property_range = self.calculate_vapor(T_range,Tc,coefs)
           path = os.path.join(tmppath,"vapor_fits",name)
        if property == "diff_cat":
           property_range = self.calculate_diffcoef(T_range,coefs)
           path = os.path.join(tmppath,"diff_cat_fits",name)
        if property == "diff_an":
           property_range = self.calculate_diffcoef(T_range,coefs)
           path = os.path.join(tmppath,"diff_an_fits",name) 
        if property == "sld":
           [property_range,T_range] = self.calculate_reduced_sld(list(T_range),Tc,coefs,epsilon_ref)
           [property_list,T_range_actual] = self.calculate_reduced_sld(list(T_range_actual),Tc,coefs,epsilon_ref)
           for i in range(len(noise_list)):
               noise_list[i] = self.reduce_density(noise_list[i],sigma_ref)
           path = os.path.join(tmppath,"guggenheim_fits",name)
           if os.getenv("METHOD") in ["stoll","stoll2"]:
              property_notred_range = []
              for i in range(len(property_range)):
                  property_notred_range.append(self.deduce_density(property_range[i],sigma_ref))
              property_notred_list = []
              for i in range(len(property_list)):
                  property_notred_list.append(self.deduce_density(property_list[i],sigma_ref))
              self.opt.append_parameter(property_notred_range,os.path.join(path,"sld_range.txt"))
              self.opt.write_vector(property_notred_list,os.path.join(path,"sld_simulated_range.txt"))
        if property == "pressure":
           property_range = self.calculate_pressure(T_range,coefs)
           path = os.path.join(tmppath,"pressure_fits",name)
        
        try:
		os.stat(path)
        except:
		os.makedirs(path)

        T_range_file = os.path.join(path,"plot_temperature.txt")
        fT = open(T_range_file,"w")
        for i in range(len(T_range)):
            if i != len(T_range)-1:
               fT.write(repr(T_range[i])+" ")
            else:
               fT.write(repr(T_range[i])+"\n")
        fT.close()

        T_file = os.path.join(path,"plot_simulated_temperature.txt")
        fT = open(T_file,"w")
        for i in range(len(T_range_actual)):
            if i != len(T_range_actual)-1:
               fT.write(repr(T_range_actual[i])+" ")
            else:
               fT.write(repr(T_range_actual[i])+"\n")
        fT.close()

        property_file = os.path.join(path,"plot_%s.txt" %(property))
        fp = open(property_file,"w")
        for i in range(len(property_range)):
            if i != len(property_range)-1:
               fp.write(repr(property_range[i])+" ")
            else:
               fp.write(repr(property_range[i])+"\n")
 
        if os.getenv("METHOD") in ["stoll","stoll2"] and property != "sld":
           self.opt.append_parameter(property_range,os.path.join(path,"%s_range.txt" %(property)))
        fp.close()
        

        property_sim_file = os.path.join(path,"plot_simulated_%s.txt" %(property))
        fs = open(property_sim_file,"w")
        for i in range(len(property_list)):
            if i != len(property_list)-1:
               fs.write(repr(property_list[i])+" ")
            else:
               fs.write(repr(property_list[i])+"\n")
        fs.close()

        property_orig_sim_file = os.path.join(path,"plot_orig_simulated_%s.txt" %(property))
        fos = open(property_orig_sim_file,"w")
        for i in range(len(old_property_list)):
            if i != len(old_property_list)-1:
               fos.write(repr(old_property_list[i])+" ")
            else:
               fos.write(repr(old_property_list[i])+"\n")
        fos.close()

        noise_sim_file = os.path.join(path,"plot_noise_%s.txt" %(property))
        fn = open(noise_sim_file,"w")
        for i in range(len(noise_list)):
            if i != len(noise_list)-1:
               fn.write(repr(noise_list[i])+" ")
            else:
               fn.write(repr(noise_list[i])+"\n")
        fn.close()

        os.system("R --slave < %s --args %s %s" %(plot_fits_script,path,property))

        return

    def diffcoef_fit(self,T_range,DC_list,dDC_list):
        """ nl regression coefficients for diffusion coefficient fits """

        print "Diffusion coefficient fit..."
        
        try:
		os.stat(path)
        except:
		os.makedirs(path)

        w = []
       
        for i in range(len(DC_list)):
            w.append(1.0/dDC_list[i])
            #w.append(1.0)

        A = []

        for T in T_range:
            A.append(-1.0/T)

        b = []
 
        for DC in DC_list:
            b.append(math.log(DC))

        f1 = open(os.path.join(path,"A.txt"),"w")
        for i in range(len(A)):
            f1.write(repr(A[i])+"\n")
        f1.close()
            
        f2 = open(os.path.join(path,"w.txt"),"w")
        for i in range(len(w)):
            f2.write(repr(w[i])+"\n")
        f2.close()    
 
        f3 = open(os.path.join(path,"b.txt"),"w")
        for i in range(len(b)):
            f3.write(repr(b[i])+"\n")
        f3.close()  

        del f1
        del f2
        del f3

        # fit by weighted linear regression
        fit_script = os.getenv("FITSCRIPT")

        os.system("R --slave < %s --args %s" %(fit_script,path))

        DC_coefs = self.opt.read_last_parameter(os.path.join(path,"coefs.txt"))

        os.system("rm %s %s %s" %(os.path.join(path,"A.txt"),os.path.join(path,"w.txt"),os.path.join(path,"b.txt")))

        return DC_coefs

    def calculate_diffcoef(self,T_range,DC_coefs):
        """ computes fitted diffusion coefficients """

        DC_list = []

        A = math.exp(DC_coefs[0])
        B = DC_coefs[1]

        for T in T_range:
            DC_list.append(A*math.exp(-B/T))

        return DC_list
 
    def fits(self,properties_list,noise_list,rho_list,drho_list,name,opt,T_range_actual=None,P_range_actual=None):
        """ fit functions for individual properties """

        tmppath = os.getenv("TMPPATH")
        properties = eval(os.getenv("PROPERTIES"))

        eff = os.getenv("efficiency")

        if eff == "y":
	   loop = int(os.getenv("loop"))-1
	   epsilon_ref_file = os.path.join(os.getenv("TMPPATH"),"epsilon_ref.%d" %(loop))
	   sigma_ref_file = os.path.join(os.getenv("TMPPATH"),"sigma_ref.%d" %(loop))

	   self.opt.epsilon_ref = self.opt.read_last_parameter(epsilon_ref_file)[0]
	   self.opt.sigma_ref = self.opt.read_last_parameter(sigma_ref_file)[0]

        if eff == None or (eff == "y" and (name.startswith("original") or name.startswith("armijo"))):
           T_range = eval(os.getenv("TRANGE"))
        else:
           T_range = eval(os.getenv("TRANGE_ERROR"))

        for prop in properties:
            if prop == "density" and rho_list == []:
               rho_index = properties.index("density")
               rho_list = []
               drho_list = []
               for k in range(len(T_range)):
                   rho_list.append(properties_list[k][rho_index])
                   drho_list.append(noise_list[k][rho_index])
            if prop == "vapor":
               vapor_index = properties.index("vapor")
               vapor_list = []
               dvapor_list = []
               for k in range(len(T_range)):
                   vapor_list.append(properties_list[k][vapor_index])
                   dvapor_list.append(noise_list[k][vapor_index])
            if prop == "diff_cat":
               DC_cat_index = properties.index("diff_cat")
               DC_cat_list = []
               dDC_cat_list = []
               for k in range(len(T_range)):
                   DC_cat_list.append(properties_list[k][DC_cat_index])
                   dDC_cat_list.append(noise_list[k][DC_cat_index])
            if prop == "diff_an":
               DC_an_index = properties.index("diff_an")
               DC_an_list = []
               dDC_an_list = []
               for k in range(len(T_range)):
                   DC_an_list.append(properties_list[k][DC_an_index])
                   dDC_an_list.append(noise_list[k][DC_an_index])
            if prop == "sld":
               sld_index = properties.index("sld")
               sld_list = []
               dsld_list = []
               for k in range(len(T_range)):
                   sld_list.append(properties_list[k][sld_index])
                   dsld_list.append(noise_list[k][sld_index])
            if prop == "pressure":
               p_index = properties.index("pressure")
               p_list = []
               dp_list = []
               for k in range(len(T_range)):
                   p_list.append(properties_list[k][p_index])
                   dp_list.append(noise_list[k][p_index])


        T_range_err = eval(os.getenv("TRANGE_ERROR"))
        if T_range_actual == None:
           T_range_actual = list(T_range_err)

        # Fits fuer saemtliche Eigenschaften
        fitted_properties = []
        if "density" in properties and os.getenv("fits") != None:
           if os.getenv("density_fit") == "specific":
              P_range = eval(os.getenv("PRANGE"))
              if P_range_actual == None:
                 P_range_actual = list(P_range)
              rho_coefs = self.density_fit(T_range_actual,P_range_actual,rho_list,drho_list,name)
           else:
              rho_coefs = self.linear_density_fit(T_range_actual,P_range_actual,rho_list,drho_list,name)
           old_rho_list = list(rho_list)
           if os.getenv("density_fit") == "specific":
              rho_list = self.calculate_density(T_range_err,rho_coefs)
              new_drho_list = []
              new_old_rho_list = []
              P0 = float(os.getenv("PRESSURE"))
              P0_index = P_range.index(P0)
              for k in range(len(T_range_err)):
                  new_drho_list.append(drho_list[P0_index][k])
                  new_old_rho_list.append(old_rho_list[P0_index][k])
              drho_list = list(new_drho_list)
              old_rho_list = list(new_old_rho_list)
              rho_list_for_reduction = []
              for P in P_range:
                  rho_list_for_reduction.append(self.calculate_density(T_range_err,rho_coefs,P))
           else:
              rho_list = self.calculate_linear_density(T_range_err,rho_coefs)
           fitted_properties.append("density")
        if "pressure" in properties and os.getenv("fits") != None:
           p_coefs = self.pressure_fit(T_range_actual,p_list,dp_list,name)
           old_p_list = list(p_list)
           p_list = self.calculate_pressure(T_range_err,p_coefs)
           fitted_properties.append("pressure")
        if "vapor" in properties and os.getenv("fits") != None:
           [vapor_coefs,Tc_vap] = self.vapor_fit(T_range_actual,vapor_list,dvapor_list,name)
           old_vapor_list = list(vapor_list)
           vapor_list = self.calculate_vapor(T_range_err,Tc_vap,vapor_coefs)
           fitted_properties.append("vapor")
        if ("diff_cat" in properties or "diff_an" in properties) and os.getenv("fits") != None:
           DC_cat_coefs = self.diffcoef_fit(T_range_actual,DC_cat_list,dDC_cat_list,name,"diff_cat")
           old_DC_cat_list = list(DC_cat_list)
           DC_cat_list = self.calculate_diffcoef(T_range_err,DC_cat_coefs)
           DC_an_coefs = self.diffcoef_fit(T_range_actual,DC_an_list,dDC_an_list,name,"diff_an")
           old_DC_an_list = list(DC_an_list)
           DC_an_list = self.calculate_diffcoef(T_range_err,DC_an_coefs)
           fitted_properties.append("diffcoef")
        if "sld" in properties and os.getenv("fits") != None:
           try:
		self.opt.sigma_ref
	   except:
		self.opt.sigma_ref = 1.0
   	   try:
		self.opt.epsilon_ref
	   except:
		self.opt.epsilon_ref = 1.0
           [sld_coefs,Tc_sld] = self.guggenheim_fit(list(T_range_actual),sld_list,dsld_list,name,self.opt.epsilon_ref,self.opt.sigma_ref)
           old_sld_list = list(sld_list)
           sld_list = self.calculate_sld(T_range_err,Tc_sld,sld_coefs,self.opt.epsilon_ref,self.opt.sigma_ref)
           fitted_properties.append("sld")

        if fitted_properties != []:
           T_range_for_plotting = self.calc_T_range_for_plotting(T_range_actual)


        for prop in fitted_properties:
            if prop == "density":
               self.plot_fits(rho_coefs,T_range_for_plotting,T_range_actual,rho_list,drho_list,name,old_rho_list,"density")
            if prop == "vapor":
               self.plot_fits(vapor_coefs,T_range_for_plotting,T_range_actual,vapor_list,dvapor_list,name,old_vapor_list,"vapor",Tc_vap)
            if prop == "diffcoef":
               self.plot_fits(DC_cat_coefs,T_range_for_plotting,T_range_actual,DC_cat_list,dDC_cat_list,name,old_DC_cat_list,"diff_cat")
               self.plot_fits(DC_an_coefs,T_range_for_plotting,T_range_actual,DC_an_list,dDC_an_list,name,old_DC_an_list,"diff_an")
            if prop == "sld":
               self.plot_fits(sld_coefs,T_range_for_plotting,T_range_actual,sld_list,dsld_list,name,old_sld_list,"sld",Tc_sld,self.opt.epsilon_ref,self.opt.sigma_ref)
            if prop == "pressure":
               self.plot_fits(p_coefs,T_range_for_plotting,T_range_actual,p_list,dp_list,name,old_p_list,"pressure")
       

        outpath = os.getenv("ORIG_OUTPATH")
        properties_file = os.path.join(outpath,"properties.txt")
        if os.getenv("METHOD") in ["stoll","stoll2"]:
           properties_range_file = os.path.join(outpath,"properties_range.txt")

        orig_properties_list = list(properties_list)

        properties_list = []
        for k in range(len(T_range_err)):
            properties_isothermal = []
            for prop in properties:
                properties_isothermal.append(0.0)
            properties_list.append(properties_isothermal)

        if os.getenv("METHOD") in ["stoll","stoll2"]:
           properties_range = []
           for k in range(len(T_range_for_plotting)):
               properties_isothermal = []
               for prop in properties:
                   properties_isothermal.append(0.0)
               properties_range.append(properties_isothermal)

        if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
           reduced_file = os.path.join(outpath,"reduced_properties.txt")

           reduced_list = []
           for k in range(len(T_range_err)):
               reduced_isothermal = []
               for prop in properties:
                   reduced_isothermal.append(0.0)
               reduced_list.append(reduced_isothermal)

           if "density" in properties and os.getenv("density_fit") == "specific":
              reduced_density_file = os.path.join(outpath,"reduced_densities.txt")
              reduced_density_list = []
              for p in range(len(P_range)):
                  reduced_density_isobaric = [] 
                  for k in range(len(T_range_err)):
                      reduced_density_isobaric.append(self.reduce_density(rho_list_for_reduction[p][k],sigma_ref))
                  reduced_density_list.append(reduced_density_isobaric)

        

        for prop in properties:
            if prop == "density":
               rho_index = properties.index("density")
               for k in range(len(T_range_err)):
                   properties_list[k][rho_index] = rho_list[k]
                   if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
                      reduced_list[k][rho_index] = self.reduce_density(rho_list[k],self.opt.sigma_ref)
               if os.getenv("METHOD") in ["stoll","stoll2"]:
                  rho_range =self.opt.read_last_parameter(os.path.join(tmppath,"density_fits",name,"density_range.txt"))
                  for k in range(len(T_range_for_plotting)):
                      properties_range[k][rho_index] = rho_range[k]
            elif prop == "vapor":
               vapor_index = properties.index("vapor")
               for k in range(len(T_range_err)):
                   properties_list[k][vapor_index] = vapor_list[k]
                   if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
                      reduced_list[k][vapor_index] = self.reduce_vapor(vapor_list[k],self.opt.epsilon_ref)
               if os.getenv("METHOD") in ["stoll","stoll2"]:
                  vapor_range =self.opt.read_last_parameter(os.path.join(tmppath,"vapor_fits",name,"vapor_range.txt"))
                  for k in range(len(T_range_for_plotting)):
                      properties_range[k][vapor_index] = vapor_range[k]
            elif prop == "diff_cat":
               DC_cat_index = properties.index("diff_cat")
               for k in range(len(T_range_err)):
                   properties_list[k][DC_cat_index] = DC_cat_list[k]
                   if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
                      reduced_list[k][DC_cat_index] = self.reduce_diffcoef(DC_cat_list[k],self.opt.epsilon_ref,self.opt.sigma_ref)
               if os.getenv("METHOD") in ["stoll","stoll2"]:
                  DC_cat_range =self.opt.read_last_parameter(os.path.join(tmppath,"diff_cat_fits",name,"diff_cat_range.txt"))
                  for k in range(len(T_range_for_plotting)):
                      properties_range[k][DC_cat_index] = DC_cat_range[k]
            elif prop == "diff_an":
               DC_an_index = properties.index("diff_an")
               for k in range(len(T_range_err)):
                   properties_list[k][DC_an_index] = DC_an_list[k]
                   if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
                      reduced_list[k][DC_an_index] = self.reduce_diffcoef(DC_an_list[k],self.opt.epsilon_ref,self.opt.sigma_ref)
               if os.getenv("METHOD") in ["stoll","stoll2"]:
                  DC_an_range =self.opt.read_last_parameter(os.path.join(tmppath,"diff_an_fits",name,"diff_an_range.txt"))
                  for k in range(len(T_range_for_plotting)):
                      properties_range[k][DC_an_index] = DC_an_range[k]
            elif prop == "sld":
               sld_index = properties.index("sld")
               for k in range(len(T_range_err)):
                   properties_list[k][sld_index] = sld_list[k]
                   if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
                      reduced_list[k][sld_index] = self.reduce_density(sld_list[k],self.opt.sigma_ref)
               if os.getenv("METHOD") in ["stoll","stoll2"]:
                  sld_range =self.opt.read_last_parameter(os.path.join(tmppath,"guggenheim_fits",name,"sld_range.txt"))
                  for k in range(len(T_range_for_plotting)):
                      properties_range[k][sld_index] = sld_range[k]
            elif prop == "pressure":
               p_index = properties.index("pressure")
               for k in range(len(T_range_err)):
                   properties_list[k][p_index] = p_list[k]
               if os.getenv("METHOD") in ["stoll","stoll2"]:
                  pressure_range =self.opt.read_last_parameter(os.path.join(tmppath,"pressure_fits",name,"pressure_range.txt"))
                  for k in range(len(T_range_for_plotting)):
                      properties_range[k][p_index] = pressure_range[k]
            else:
               prop_index = properties.index(prop)
               for k in range(len(T_range_err)):
                   properties_list[k][prop_index] = orig_properties_list[k][prop_index]


        f = open(properties_file,"w")
        for k in range(len(T_range_err)):
            for l in range(len(properties)):
                f.write(repr(properties_list[k][l])+"\n")
        f.close()

        if os.getenv("METHOD") in ["stoll","stoll2"]:
           f = open(properties_range_file,"w")
           for k in range(len(T_range_for_plotting)):
               for l in range(len(properties)):
                   f.write(repr(properties_range[k][l])+"\n")
           f.close()

        if eff == "y" and (name.startswith("original") or name.startswith("armijo")):
           f_red = open(reduced_file,"w")
           for k in range(len(T_range_err)):
               for l in range(len(properties)):
                   if l != len(properties)-1:
                      f_red.write(repr(reduced_list[k][l])+" ")
                   else:
                      f_red.write(repr(reduced_list[k][l])+"\n")
           f_red.close()
  
           if "density" in properties and os.getenv("density_fit") == "specific":
              f_dens_red = open(reduced_density_file,"w")
              for p in range(len(P_range)):
                  for k in range(len(T_range_err)):
                      if k != len(T_range_err)-1:
                         f_dens_red.write(repr(reduced_density_list[p][k])+" ")
                      else:
                         f_dens_red.write(repr(reduced_density_list[p][k])+"\n") 
              f_dens_red.close()

        return


    def reduce_parameters(self,x,dim):
        """ calculates the reduced parameters """

        par_red_list = []

        param_order = eval(os.getenv("PARAMORDER"))
        param_number = eval(os.getenv("PARAMNUMBER"))
        parind_list = eval(os.getenv("PARINDLIST"))

        if "epsilon" in param_order:
           eps_ind = param_order.index("epsilon")
           if eps_ind == 0:
              eps_start = 0
           elif eps_ind == 1:
              eps_start = param_number[0]
           else:
              eps_start = param_number[0]+param_number[1]
           eps_ref = x[eps_start]
        else:
           eps_ref = None

        if "sigma" in param_order:
           sigma_ind = param_order.index("sigma")
           if sigma_ind == 0:
              sigma_start = 0
           elif sigma_ind == 1:
              sigma_start = param_number[0]
           else:
              sigma_start = param_number[0]+param_number[1]
           sigma_ref = x[sigma_start]
        else:
           sigma_ref = None

        if "charge" in param_order:
           charge_ind = param_order.index("charge")
           if charge_ind == 0:
              charge_start = 0
           elif charge_ind == 1:
              charge_start = param_number[0]
           else:
              charge_start = param_number[0]+param_number[1]
           charge_red_factor = float(os.getenv("charge_red_factor"))
           eps_ref = x[eps_start]
           charge_factor = charge_red_factor/eps_ref

        for i in range(dim):
            parameter =self.opt.which_parameter(param_order,parind_list,i+1)
            if parameter == "epsilon":
               if i == eps_start:
                  par_red_list.append(1.0)
               else:
                  par_red_list.append(x[i]/eps_ref)
            if parameter == "sigma":
               if i == sigma_start:
                  par_red_list.append(1.0)
               else:
                  par_red_list.append(x[i]/sigma_ref)
            if parameter == "charge":              
               par_red_list.append(charge_factor*x[i])

        return [par_red_list,eps_ref,sigma_ref]

    def change_parameter(self,y,h,i,par_red_list):
        """ changes force field parameters simultaneously """

        param_order = eval(os.getenv("PARAMORDER"))
        param_number = eval(os.getenv("PARAMNUMBER"))
        parind_list = eval(os.getenv("PARINDLIST"))

        parameter =self.opt.which_parameter(param_order,parind_list,i+1)

        y[i] = y[i] + h

        if parameter == "epsilon" or parameter == "sigma":
           par_ind = param_order.index(parameter)
           if par_ind == 0:
              par_start = 0
           elif par_ind == 1:
              par_start = param_number[0]
           else:
              par_start = param_number[0]+param_number[1]
           nof_pars = param_number[par_ind]
           
           # andere epsilons/sigmas aendern sich mit
           indic = i-par_start
           par_ref_tilde = y[i]/par_red_list[i] # beachte: y[i] = y[i] + h
           for j in range(nof_pars):
               if j != indic:
                  y[par_start+j] = par_ref_tilde*par_red_list[par_start+j]

        if parameter == "epsilon" and "charge" in param_order:
           # Ladungen veraendern sich bei Aenderung von epsilon !!!
           charge_ind = param_order.index("charge")
           if charge_ind == 0:
              charge_start = 0
           elif charge_ind == 1:
              charge_start = param_number[0]
           else:
              charge_start = param_number[0]+param_number[1]
           nof_charges = param_number[charge_ind]
            
           par_ref = y[par_start]/par_red_list[par_start]
 
           for j in range(nof_charges):
               y[charge_start+j] = par_ref_tilde/par_ref * y[charge_start+j]
           
           

        if parameter == "charge":
           # epsilons veraendern sich bei Aenderung der Ladung und somit auch alle anderen Ladungen !!!
           charge_red_factor = float(os.getenv("charge_red_factor"))
           charge_ind = param_order.index("charge")
           if charge_ind == 0:
              charge_start = 0
           elif charge_ind == 1:
              charge_start = param_number[0]
           else:
              charge_start = param_number[0]+param_number[1]
           nof_charges = param_number[charge_ind]
           
           eps_ind = param_order.index("epsilon")
           if eps_ind == 0:
              eps_start = 0
           elif eps_ind == 1:
              eps_start = param_number[0]
           else:
              eps_start = param_number[0]+param_number[1]
           nof_eps = param_number[eps_ind]

           par_ref_tilde = y[i]/par_red_list[i] * charge_red_factor

           # Aenderung der epsilons
           for j in range(nof_eps):
               y[eps_start+j] = par_ref_tilde*par_red_list[eps_start+j]

           par_ref = y[eps_start]/par_red_list[eps_start]

           # Aenderung der anderen Ladungen
           indic = i-charge_start
           for j in range(nof_charges):
               if j != indic:
                  y[charge_start+j] = par_ref_tilde/par_ref * y[charge_start+j]
                

        self.opt.append_parameter(y,os.getenv("PARAMFILE"))

        return [parameter,par_ref_tilde]

    def change_two_parameters(self,y,h,i,j,par_red_list):
        """ changes force field parameters simultaneously for Hessian matrix """

        param_order = eval(os.getenv("PARAMORDER"))
        param_number = eval(os.getenv("PARAMNUMBER"))
        parind_list = eval(os.getenv("PARINDLIST"))

        parameter1 =self.opt.which_parameter(param_order,parind_list,i+1)
        parameter2 =self.opt.which_parameter(param_order,parind_list,j+1)

        y[i] = y[i] + h
        y[j] = y[j] + h

        sigma_eps = False
        sigma_sigma = False
        eps_eps = False
       
        sigma_charge = False
        eps_charge = False
        charge_charge = False

        if parameter1 == "sigma" and parameter2 == "epsilon":
           sigma_ind = param_order.index(parameter1)
           eps_ind = param_order.index(parameter2)
           sigma_eps = True
        if parameter1 == "epsilon" and parameter2 == "sigma":
           sigma_ind = param_order.index(parameter2)
           eps_ind = param_order.index(parameter1)
           sigma_eps = True
        if parameter1 == "sigma" and parameter2 == "sigma":
           sigma_ind = param_order.index(parameter1)
           sigma_sigma = True
        if parameter1 == "epsilon" and parameter2 == "epsilon":
           eps_ind = param_order.index(parameter1)
           eps_eps = True
        if parameter1 == "sigma" and parameter2 == "charge":
           sigma_ind = param_order.index(parameter1)
           charge_ind = param_order.index(parameter2)
           sigma_charge = True
        if parameter1 == "charge" and parameter2 == "sigma": 
           sigma_ind = param_order.index(parameter2)
           charge_ind = param_order.index(parameter1)
           sigma_charge = True
        if parameter1 == "epsilon" and parameter2 == "charge":
           eps_ind = param_order.index(parameter1)
           charge_ind = param_order.index(parameter2)
           eps_charge = True
        if parameter1 == "charge" and parameter2 == "epsilon": 
           eps_ind = param_order.index(parameter2)
           charge_ind = param_order.index(parameter1)
           eps_charge = True
        if parameter1 == "charge" and parameter2 == "charge":
           charge_ind = param_order.index(parameter1)
           charge_charge = True

        if sigma_eps == True or sigma_sigma == True or sigma_charge == True:
           if sigma_ind == 0:
              sigma_start = 0
           elif sigma_ind == 1:
              sigma_start = param_number[0]
           else:
              sigma_start = param_number[0]+param_number[1]
           nof_sigmas = param_number[sigma_ind]
           

        if sigma_eps == True or eps_eps == True or eps_charge == True or sigma_charge == True or charge_charge == True:
           if eps_ind == 0:
              eps_start = 0
           elif eps_ind == 1:
              eps_start = param_number[0]
           else:
              eps_start = param_number[0]+param_number[1]          
           nof_epsilons = param_number[eps_ind]
          
        if (sigma_charge == True or eps_charge == True or charge_charge == True) or ((sigma_eps == True or eps_eps == True) and "charge" in param_order):
           charge_red_factor = float(os.getenv("charge_red_factor"))
           charge_ind = param_order.index("charge")
           if charge_ind == 0:
              charge_start = 0
           elif charge_ind == 1:
              charge_start = param_number[0]
           else:
              charge_start = param_number[0]+param_number[1]
           nof_charges = param_number[charge_ind]

        if sigma_eps == True:
           # andere epsilons/sigmas aendern sich mit
           if parameter1 == "epsilon" and parameter2 == "sigma":
              sigma_indic = j-sigma_start
              eps_indic = i-eps_start
              sigma_ref_tilde = y[j]/par_red_list[j] # beachte: y[j] = y[j] + h
              eps_ref_tilde = y[i]/par_red_list[i] # beachte: y[i] = y[i] + h
           else:
              sigma_indic = i-sigma_start
              eps_indic = j-eps_start
              sigma_ref_tilde = y[i]/par_red_list[i] # beachte: y[i] = y[i] + h
              eps_ref_tilde = y[j]/par_red_list[j] # beachte: y[j] = y[j] + h
           for k in range(nof_sigmas):
               if k != sigma_indic:
                  y[sigma_start+k] = sigma_ref_tilde*par_red_list[sigma_start+k]
           for k in range(nof_epsilons):
               if k != eps_indic:
                  y[eps_start+k] = eps_ref_tilde*par_red_list[eps_start+k]
           if parameter1 == "epsilon" and parameter2 == "sigma":
              par1_ref_tilde = eps_ref_tilde
              par2_ref_tilde = sigma_ref_tilde
           else:
              par1_ref_tilde = sigma_ref_tilde
              par2_ref_tilde = eps_ref_tilde

        if sigma_sigma == True:
           # andere sigmas aendern sich mit
           sigma1_indic = i-sigma_start
           sigma2_indic = j-sigma_start
           sigma_ref_tilde = (y[i] + par_red_list[i]/par_red_list[j])/par_red_list[i] # beachte: y[i] = y[i] + h
           for k in range(nof_sigmas):
               if k != sigma1_indic and k != sigma2_indic:
                  y[sigma_start+k] = sigma_ref_tilde*par_red_list[sigma_start+k]
           par1_ref_tilde = sigma_ref_tilde
           par2_ref_tilde = sigma_ref_tilde

        if eps_eps == True:
           # andere epsilons aendern sich mit
           eps1_indic = i-eps_start
           eps2_indic = j-eps_start
           eps_ref_tilde = (y[i] + par_red_list[i]/par_red_list[j])/par_red_list[i] # beachte: y[i] = y[i] + h
           for k in range(nof_epsilons):
               if k != eps1_indic and k != eps2_indic:
                  y[eps_start+k] = eps_ref_tilde*par_red_list[eps_start+k]
           par1_ref_tilde = epsilon_ref_tilde
           par2_ref_tilde = epsilon_ref_tilde
           
          
        if (sigma_eps == True or eps_eps == True) and "charge" in "param_order":
           # Ladungen aendern sich bei Aenderung von epsilon mit               
           eps_ref = y[eps_start]/par_red_list[eps_start]
 
           for k in range(nof_charges):
               y[charge_start+k] = eps_ref_tilde/eps_ref * y[charge_start+k]
 

        if sigma_charge == True:
           # epsilons, andere sigmas und andere Ladungen aendern sich mit

           eps_ref = y[eps_start]/par_red_list[eps_start]
          
           if parameter1 == "sigma" and parameter2 == "charge":
              sigma_indic = i-sigma_start
              charge_indic = j-charge_start
              sigma_ref_tilde = y[i]/par_red_list[i] # beachte: y[i] = y[i] + h
              eps_ref_tilde = y[j]/par_red_list[j] * charge_red_factor
           else:
              sigma_indic = j-sigma_start
              charge_indic = i-charge_start
              sigma_ref_tilde = y[j]/par_red_list[j] # beachte: y[j] = y[j] + h
              eps_ref_tilde = y[i]/par_red_list[i] * charge_red_factor

           # Aenderungen der sigmas
           for k in range(nof_sigmas):
               if k != sigma_indic:
                  y[sigma_start+k] = sigma_ref_tilde*par_red_list[sigma_start+k]
           
           # Aenderung der anderen Ladungen
           for k in range(nof_charges):
               if k != charge_indic:
                  y[charge_start+k] = eps_ref_tilde/eps_ref * y[charge_start+k]

           # Aenderung der epsilons
           for k in range(nof_eps):
               y[eps_start+k] = eps_ref_tilde*par_red_list[eps_start+k]

           if parameter1 == "sigma" and parameter2 == "charge":
              par1_ref_tilde = sigma_ref_tilde
              par2_ref_tilde = eps_ref_tilde
           else:
              par1_ref_tilde = eps_ref_tilde
              par2_ref_tilde = sigma_ref_tilde

        if eps_charge == True:
           # andere epsilons und andere Ladungen aendern sich

           eps_ref = y[eps_start]/par_red_list[eps_start]

           if parameter1 == "epsilon" and parameter2 == "charge":
              epsilon_indic = i-epsilon_start
              charge_indic = j-charge_start
              eps_ref_tilde = (y[i] + par_red_list[i]/par_red_list[j])/par_red_list[i] # beachte: y[i] = y[i] + h
           else:
              epsilon_indic = j-epsilon_start
              charge_indic = i-charge_start
              eps_ref_tilde = (y[i] + par_red_list[i]/par_red_list[j])/par_red_list[i] * charge_red_factor # beachte: y[i] = y[i] + h
           
           # Aenderung der anderen Ladungen
           for k in range(nof_charges):
               if k != charge_indic:
                  y[charge_start+k] = eps_ref_tilde/eps_ref * y[charge_start+k]

           # Aenderung der epsilons
           for k in range(nof_eps):
               if k != epsilon_indic:
                  y[eps_start+k] = eps_ref_tilde*par_red_list[eps_start+k]

           par1_ref_tilde = eps_ref_tilde
           par2_ref_tilde = eps_ref_tilde

        if charge_charge == True:
           # epsilons und andere Ladungen aendern sich

           eps_ref = y[eps_start]/par_red_list[eps_start]

           charge1_indic = i-charge_start
           charge2_indic = j-charge_start
           eps_ref_tilde = (y[i] + par_red_list[i]/par_red_list[j])/par_red_list[i] * charge_red_factor # beachte: y[i] = y[i] + h
 
           # Aenderung der anderen Ladungen
           for k in range(nof_charges):
               if k != charge1_indic and k != charge2_indic:
                  y[charge_start+k] = eps_ref_tilde/eps_ref * y[charge_start+k]

           # Aenderung der epsilons
           for k in range(nof_eps):
                  y[eps_start+k] = eps_ref_tilde*par_red_list[eps_start+k]

           par1_ref_tilde = eps_ref_tilde
           par2_ref_tilde = eps_ref_tilde

        return [parameter1,parameter2,par1_ref_tilde,par2_ref_tilde]

    def reduce_density(self,rho,sigma_ref):
        """ computes rho* """

        rho_red_factor = float(os.getenv("rho_red_factor"))

        return rho_red_factor * rho * pow(sigma_ref,3)

    def deduce_density(self,rho_red,sigma_ref):
        """ computes rho from rho* """

        rho_red_factor = float(os.getenv("rho_red_factor"))

        return rho_red / (rho_red_factor*pow(sigma_ref,3))

    def reduce_vapor(self,vapor,epsilon_ref):
        """ computes Delta hv* """

        vapor_red_factor = float(os.getenv("vapor_red_factor"))

        return vapor_red_factor * vapor/epsilon_ref

    def reduce_diffcoef(self,diffcoef,epsilon_ref,sigma_ref):
        """ computes D* """

        diffcoef_red_factor = float(os.getenv("diffcoef_red_factor"))

        return diffcoef_red_factor * diffcoef/(sigma_ref*math.sqrt(epsilon_ref))

    def reduce_temperature(self,T,epsilon_ref):
        """ computes T* """

        T_red_factor = float(os.getenv("T_red_factor"))

        return T_red_factor * T/epsilon_ref

    def reduce_pressure(self,pressure,epsilon_ref,sigma_ref):
        """ computes P* """

        pressure_red_factor = float(os.getenv("pressure_red_factor"))

        return pressure_red_factor * pressure * pow(sigma_ref,3)/epsilon_ref

    def density(self,rho_red,sigma_ref):
        """ computes rho from rho* """

        rho_red_factor = float(os.getenv("rho_red_factor"))

        return rho_red/(rho_red_factor*pow(sigma_ref,3))

    def vapor(self,vapor_red,epsilon_ref):
        """ computes Delta hv from Delta hv* """

        vapor_red_factor = float(os.getenv("vapor_red_factor"))

        return vapor_red * epsilon_ref/vapor_red_factor

    def diffcoef(self,diffcoef_red,epsilon_ref,sigma_ref):
        """ computes D from D* """

        diffcoef_red_factor = float(os.getenv("diffcoef_red_factor"))

        return diffcoef_red * sigma_ref * math.sqrt(epsilon_ref)/diffcoef_red_factor

    def temperature(self,temperature_red,epsilon_ref):
        """ computes T from T* """

        T_red_factor = float(os.getenv("T_red_factor"))

        return temperature_red * epsilon_ref/T_red_factor

    def pressure(self,pressure_red,epsilon_ref,sigma_ref):
        """ computes P from P* """

        pressure_red_factor = float(os.getenv("pressure_red_factor"))

        return pressure_red/pressure_red_factor * epsilon_ref/pow(sigma_ref,3)  
        
    def actual_T_range(self,T_range,epsilon_ref,epsilon_tilde_ref):
        """ comuptes the actual temperatures of the simulated simulation """

        for k in range(len(T_range)):
            T_range[k] = epsilon_tilde_ref/epsilon_ref * T_range[k]

        return T_range

    def actual_P_range(self,P_range,epsilon_ref,sigma_ref,epsilon_tilde_ref,sigma_tilde_ref):
        """ comuptes the actual pressures of the simulated simulation """

        for p in range(len(P_range)):
            P_range[p] = (pow(sigma_ref,3)*epsilon_tilde_ref) / (pow(sigma_tilde_ref,3)*epsilon_ref) * P_range[p]

        return P_range

    def compute_properties(self,red_par_file,epsilon_ref,sigma_ref,T_range,properties):
        """ reads the reduced properties from a file and coordinates the computations of new properties """

        properties_list = self.opt.read_matrix(red_par_file)

        if "density" in properties and os.getenv("density_fit") == "linear":
           rho_index = properties.index("density")
           rho_red = []
           for k in range(len(T_range)):
               properties_list[k][rho_index] = self.density(properties_list[k][rho_index],sigma_ref)

        if "sld" in properties:
           sld_index = properties.index("sld")
           sld_red = []
           for k in range(len(T_range)):
               properties_list[k][sld_index] = self.density(properties_list[k][sld_index],sigma_ref)

        if "vapor" in properties:
           vapor_index = properties.index("vapor")
           vapor_red = []
           for k in range(len(T_range)):
               properties_list[k][vapor_index] = self.vapor(properties_list[k][vapor_index],epsilon_ref)

        if "diff_cat" in properties:
           diff_cat_index = properties.index("diff_cat")
           diff_cat_red = []
           for k in range(len(T_range)):
               properties_list[k][diff_cat_index] = self.diffcoef(properties_list[k][diff_cat_index],epsilon_ref,sigma_ref)

        if "diff_an" in properties:
           diff_an_index = properties.index("diff_an")
           diff_an_red = []
           for k in range(len(T_range)):
               properties_list[k][diff_an_index] = self.diffcoef(properties_list[k][diff_an_index],epsilon_ref,sigma_ref)
        
        noise_list = []

        for k in range(len(T_range)):
            noise_isothermal = []
            for l in range(len(properties)):
                noise_isothermal.append(1.0)
            noise_list.append(noise_isothermal)

        return [properties_list,noise_list]

    def compute_density(self,red_density_file,sigma_ref,P_range,T_range):
        """ reads the reduced densities from a file and computes the new densities """

        density_matrix = self.opt.read_matrix(red_density_file)

        for p in range(len(P_range)):
            for k in range(len(T_range)):
                density_matrix[p][k] = self.density(density_matrix[p][k],sigma_ref)

        noise_density_matrix = []

        for p in range(len(P_range)):
            noise_isobaric = []
            for k in range(len(T_range)):
                noise_isobaric.append(1.0)
            noise_density_matrix.append(noise_isobaric)

        return [density_matrix,noise_density_matrix]

    def __del__(self):
          """
          destructor
          """	
	  del self







































