#  class Physical_Properties_Loss, 2017/02/07 huelsmann FhI SCAI/HBRS
#
#start_doc
#
#  Script:          physical_properties_loss.py
#  
#  Author:          Marco Huelsmann
#                      Robin Strickstrock
# 
#  Date:             07-02-2017
#
#  Description:    class Physical_Properties_Loss (derived from Loss
#                      combining former classes:
#                      physical_properties_loss
#                    & qm_mm_loss
#  Usage:            by defining an instance
#
#  Arguments:
#
#  Options:      
# 			
#  Output:         
#
#  Imported:        loss.py
#                      opt_class.py
#
#  Calling:         parameter_substance_ensemble_temperature_variation.py     (physical_properties_loss)
#                      parameter_variation.py                                              (qm_mm_loss)
#                      executer.py
#                      collector.py                                                            (physical_properties_loss)
# 			
#
#  Modifications: hybrid functionality for physical_properties_loss & qm_mm_loss combined
#
#  $Log: md_opt$
#
#end_doc

# ***********************************************************************
# required python modules
# ***********************************************************************

from generic_optimization_problem.loss import Loss

import os
import os.path
import string
import sys
import subprocess
import time

from utilities.io import IO
i_o = IO()
from utilities.math_utilities import Math
m = Math()
from utilities.trace import Trace
tr = Trace()
from utilities.system import System
sy = System()
from utilities.string_utilities import String_Utilities
su = String_Utilities()

## only used by Physical_Properties_Loss ##
import opt_class
opt = opt_class.Optimization()


# ***********************************************************************
# class Physical_Properties_Loss
# ***********************************************************************

class Physical_Properties_Loss(Loss):
         
    def __init__(self, dimension, weights, targets, config):
        """ constructor """

        self.__config = config 
        self.__cwd = os.path.dirname(os.path.realpath(sys.argv[0]))
        #print self.__cwd

        ### checking objective_function in config-file ###
        self._objective_function = self.__config.get("OPT","objective_function")
        if self._objective_function not in ["QM_MM_Loss","Physical_Properties_Loss","PhysProp_QMMM_Loss"]:
            print "ERROR in physical_properties_loss_new.__init__: invalid objective_function (%s)! Check config file." %(self._objective_function)
            tr.errorexit()        

        sections = list(self.__config.sections())

        if self._objective_function == "QM_MM_Loss" or self._objective_function == "Physical_Properties_Loss":
            self._targets = []
            self._weights = []

            for section in sections:
                if section not in ["SYS", "OPT", "META"]:
                    self.__init_targets(section)

            for section in sections:
                if section not in ["SYS", "OPT", "META"]:
                    self.__init_weights(section)

        if self._objective_function == "PhysProp_QMMM_Loss":
            self._targets_QM_MM = []
            self._targets_Physical_Properties = []

            self._weights_QM_MM = []
            self._weights_Physical_Properties = []

            for section in sections:
                if section not in ["SYS", "OPT", "META"]:
                    self.__init_targets(section)
            self._targets = self._targets_QM_MM
            for target_PP in self._targets_Physical_Properties:
                self._targets.append(target_PP)

            for section in sections:
                if section not in ["SYS", "OPT", "META"]:
                    self.__init_weights(section)
            self._weights = self._weights_QM_MM
            for weight_PP in self._weights_Physical_Properties:
                self._weights.append(weight_PP)

        # weight_mode
        try:
            weight_mode = self.__config.get("OPT", "weight_mode")
        except:
            print "You must indicate how to define the weights."
            tr.errorexit()

        if weight_mode not in ["as_is", "compute", "normalize"]:
            print "For 'weight_mode' you must choose between 'as_is' and 'compute'. %s is not valid." % (weight_mode)
            tr.errorexit()

        if weight_mode == "compute" or weight_mode == "normalize":
            weight_sum = sum(self._weights)
            #print "weight_sum = %s" %(str(weight_sum))
            for i in range(len(self._weights)):
                #print "i = %s" %(str(i))
                #print "weight[%s] = %s" %(str(i), self._weights[i])
                self._weights[i] /= weight_sum
                #print "new weight[%s] = %s" %(str(i),self._weights[i]) 
            #print "new weight_sum = %s" %(str(sum(self._weights)))

        Loss.__init__(self, dimension, self._weights, self._targets)


    def __init_targets(self, section):
        """ initializes the targets """

        # target          
        try:
            os.stat(self.__config.get(section, "target"))
        except OSError:
            print "ERROR in physical_properties_loss_new.__init_targets: Indicated target file %s does not exist." % (self.__config.get(section, "target"))
            tr.errorexit()
        except:
            print "ERROR in physical_properties_loss_new.__init_targets: You must indicate a target file."
            tr.errorexit()
        else:
            if self._objective_function == "PhysProp_QMMM_Loss":
                # section for physical properties loss
                if section in ["MD"]:
                    target_file_Physical_Properties = os.path.abspath(self.__config.get(section, "target"))
                    #print "Physical_Properties_Loss_new: target_file_Physical_Properties:",target_file_Physical_Properties
                # section for qm_mm loss
                if section in ["QM_MM"]:
                    target_file_QM_MM = os.path.abspath(self.__config.get(section, "target"))
                    #print "Physical_Properties_Loss_new: target_file_QM_MM:",target_file_QM_MM
                if section not in ["MD","QM_MM"]:
                    print "ERROR in pyhsical_properties_loss_new.__init_targets: Only use sections [MD] and [QM_MM] (besides [SYS], [OPT] and [META]). Section [%s] is not supported." %(section)
                    tr.errorexit()  

            if self._objective_function == "QM_MM_Loss" or self._objective_function == "Physical_Properties_Loss":
                target_file = os.path.abspath(self.__config.get(section, "target"))

        if self._objective_function == "QM_MM_Loss" or self._objective_function == "Physical_Properties_Loss":
            target_file = os.path.abspath(target_file)
            targets = i_o.get_properties(target_file)
            for target in targets:
                self._targets.append(target)

            sy.setenv("TARGET", target_file)

        if self._objective_function == "PhysProp_QMMM_Loss":
            # section for physical properties loss
            if section in ["MD"]:
                target_file_Physical_Properties = os.path.abspath(target_file_Physical_Properties)
                targets_Physical_Properties = i_o.get_properties(target_file_Physical_Properties)
                for target_Physical_Properties in targets_Physical_Properties:
                    self._targets_Physical_Properties.append(target_Physical_Properties)

                sy.setenv("TARGET_Physical_Properties", target_file_Physical_Properties)

            # section for qm_mm loss
            if section in ["QM_MM"]:
                target_file_QM_MM = os.path.abspath(target_file_QM_MM)
                targets_QM_MM = i_o.get_properties(target_file_QM_MM)
                for target_QM_MM in targets_QM_MM:
                    self._targets_QM_MM.append(target_QM_MM)

                sy.setenv("TARGET_QM_MM", target_file_QM_MM)


    def __init_weights(self, section):
        """ initializes the weights """        

        ### QM_MM_Loss ###
        if self._objective_function == "QM_MM_Loss":
            self.__init_weights_QM_MM()

        ### Physical_Properties_Loss ###
        if self._objective_function == "Physical_Properties_Loss":
            self.__init_weights_Physical_Properties(section)

        ### Physical_Properties_Loss and QM_MM_Loss combined ###
        if self._objective_function == "PhysProp_QMMM_Loss":

            if section in ["MD"]:
                self.__init_weights_Physical_Properties(section)
            if section in ["QM_MM"]:
                self.__init_weights_QM_MM()


    def __init_weights_QM_MM(self):
        weight = float(self.__config.get("QM_MM","weights"))
        #print "in physical_properties_loss_new.py: weight=",weight
        try:
            weight_scaling = self.__config.get("QM_MM", "weight_scaling")
        except:
            print "weight_scaling in %s [QM_MM] not set. Every QM_MM target is weighted 'as_is' or 'compute'."
            weight_scaling = None


        if self._objective_function == "QM_MM_Loss":
            for i in range(len(self._targets)):
                self._weights.append(weight)

        if self._objective_function == "PhysProp_QMMM_Loss":
            i=1
            for i in range(len(self._targets_QM_MM)-1):
                if weight_scaling == "ignore_lower_half":
                    if i < float(len(self._targets_QM_MM))/2:
                        self._weights_QM_MM.append(0)
                    else:
                        self._weights_QM_MM.append(weight*2)
                    i+=1

                elif weight_scaling == "ignore_higher_half":
                    if i < float(len(self._targets_QM_MM))/2:
                        self._weights_QM_MM.append(weight*2)
                    else:
                        self._weights_QM_MM.append(0)
                    i+=1
                else:
                    self._weights_QM_MM.append(weight)
            #print "in physical_properties_loss_new.py: self._weights_QM_MM:", self._weights_QM_MM


    def __init_weights_Physical_Properties(self, section):
        # properties          
        try:  
            properties = self.__config.get(section, "properties")
        except:
            print "ERROR in physical_properties_loss_new.__init_weights: You must indicate the physical properties considered in the optimization."
            tr.errorexit()

        properties = string.split(properties)

        # weights
        try:
            weights = self.__config.get(section, "weights")
        except:
            print "ERROR in physical_properties_loss_new.__init_weights: You must indicate the weights for the loss function."
            tr.errorexit()

        weights = string.split(weights)
        weights = opt.eval_vector(weights, check=True)
        if type(weights) != list:
            print "ERROR in physical_properties_loss_new.__init_weights: The indicated weight '%s' is not evaluable." % (weights)
            tr.errorexit()

        # temp
        VLE_T_range = []
        T_range_trans = []
        NPT_T_range = []
        NVT_T_range = []

        NPT_VLE = "n"
        NPT_Trans = "n"
        NVT_Trans = "n"

        temp_counter = 0

        for prop in properties:

            weight_index = properties.index(prop)

    ### not tested
    ####### temp for 'sld' | 'dhv' | 'ps' ##########
            if prop == "sld" or prop == "dhv" or prop == "ps":
                if VLE_T_range == []:
                    try:
                        VLE_T = self.__config.get(section, "VLE_temperature")
                    except:
                        print "You must indicate the VLE temperature range for the loss function."
                        tr.errorexit()

                    VLE_T = string.split(VLE_T)
                    VLE_T_range = su.float_vector(VLE_T, check=True)
                    if type(VLE_T_range) != list:
                        print "The indicated VLE temperature '%s' is not numeric." % (VLE_T_range)
                        tr.errorexit()

                    temp = list(VLE_T_range)
    ### not tested
    ####### temp for 'diff_cat' | 'diff_an' | 'reorient_time' | 'viscosity' | 'thermal_cond' ##########
            if prop == "diff_cat" or prop == "diff_an" or prop == "reorient_time" or prop == "viscosity" or prop == "thermal_cond":
                if T_range_trans == []:
                    try:
                        T_range_trans = self.__config.get(section, "Trans_temperature")
                    except:
                        print "You must indicate the temperature range ('Trans_temperature') for the calculation of transport properties."
                        tr.errorexit()

                    T_range_trans = string.split(T_range_trans)
                    T_range_trans = su.float_vector(T_range_trans, check=True)

                temp = list(T_range_trans)

    ####### temp for 'density' | 'isothermal_compress' | 'volume_expans' | 'dHdP' | 'speed' | 'isobaric_heatcap' ##########
            if prop == "density" or prop == "isothermal_compress" or prop == "volume_expans" or prop == "dHdP" or prop == "speed" or prop == "isobaric_heatcap":
                if NPT_T_range == []:
                    try:
                        NPT_T = self.__config.get(section, "NPT_temperature")
                    except:
                        if NPT_VLE == "n" and NPT_Trans == "n":
                            print "WARNING: You have indicated properties for which an NPT simulation is required and you did not indicate specific temperatures."

                        if NPT_VLE == "n":
                            print "Taking VLE temperatures NPT simulations."

                        if len(VLE_T_range) == 0:                                    
                            try:
                                VLE_T = self.__config.get(section, "VLE_temperature")
                            except:
                                print "VLE temperature list is also empty. Taking VLE temperatures."
                                try:
                                    T_range_trans = self.__config.get(section, "Trans_temperature")
                                except:
                                    print "VLE temperature list is also empty. Taking VLE temperatures."
                                    try:
                                        T_range_trans = self.__config.get(mol, "Trans_temperature")
                                    except:
                                        print "Transport temperature list is also empty. Abort now."
                                        tr.errorexit()
                                    else:
                                        T_range_trans = string.split(T_range_trans)
                                        T_range_trans = su.float_vector(T_range_trans, check=True)
                                else:                                     
                                    T_range_trans = string.split(T_range_trans)
                                    T_range_trans = su.float_vector(T_range_trans, check=True)
                            else:
                                VLE_T = string.split(VLE_T)
                                VLE_T_range = su.float_vector(VLE_T, check=True)
                                if type(VLE_T_range) != list:
                                    print "The indicated VLE temperature '%s' is not numeric." % (VLE_T_range)
                                    tr.errorexit()

                        if len(VLE_T_range) > 0:
                            NPT_T_range = list(VLE_T_range)
                            NPT_VLE = "y"

                        if len(T_range_trans) > 0 and len(NPT_T_range) == 0:
                            NPT_T_range = list(T_range_trans)
                            NPT_Trans = "y"
                    else:
                        NPT_T = string.split(NPT_T)
                        NPT_T_range = su.float_vector(NPT_T, check=True)
                        if type(NPT_T_range) != list:
                            print "The indicated NPT temperature '%s' is not numeric." % (NPT_T_range)
                            tr.errorexit()

                temp = list(NPT_T_range)

    ### not tested
    ####### temp for 'pressure' | 'isochoric_heatcap' | 'dUdV' ##########
            if prop == "pressure" or prop == "isochoric_heatcap" or prop == "dUdV":
                if NVT_T_range == []:
                    try:
                        NVT_T = self.__config.get(section, "NPT_temperature") 
                    except:
                        if NVT_Trans == "n":
                            print "WARNING: You have indicated properties for which an NVT simulation is required and you did not indicate specific temperatures."
                            print "Taking transport temperatures and pressures for transport simulations."

                        if len(T_range_trans) == 0:
                            try:
                                T_range_trans = self.__config.get(section, "Trans_temperature")
                            except:
                                print "Transport temperature list is also empty. Abort now."
                                tr.errorexit()
                            else:
                                T_range_trans = string.split(T_range_trans)
                                T_range_trans = su.float_vector(T_range_trans, check=True)

                        NVT_T_range = list(T_range_trans)
                        NVT_Trans = "y"			
                    else:
                        NVT_T = string.split(NVT_T)
                        NVT_T_range = su.float_vector(NVT_T, check=True)
                        if type(NVT_T_range) != list:
                            print "The indicated NVT temperature '%s' is not numeric." % (NVT_T_range)
                            tr.errorexit()

                temp = list(NVT_T_range)
    #######
            temp_counter += len(temp)
 
            # collocate property-specific weights
            if len(weights) == len(properties):      
                if self._objective_function == "Physical_Properties_Loss":        
                    for T in temp:
                        self._weights.append(weights[weight_index])
                if self._objective_function == "PhysProp_QMMM_Loss":
                    for T in temp:
                        self._weights_Physical_Properties.append(weights[weight_index])

        # collocate property-specific weights
        if len(weights) != len(properties):
            if len(weights) != temp_counter:
                print "ERROR in physical_properties_loss_new.__init_weights: The number of weights (%d) does not correspond either to the number of properties (%d) nor to the number of temperatures (%d)!" %(len(weights), len(properties),temp_counter)
                tr.errorexit()
            else:
                if self._objective_function == "Physical_Properties_Loss":
                    for w in weights:
                        self._weights.append(w)
                if self._objective_function == "PhysProp_QMMM_Loss":
                    for w in weights:
                        self._weights_Physical_Properties.append(w)


    def get_function_values(self, parameter_set, begin, return_estimations = False):
        """ returns a list of function values for a given parameter set"""
          
        sy.setenv("BEGIN", repr(begin))

        parameter_set_file = os.path.join(os.getenv("GROW_HOME"), "parameter_set.%s" %(os.getenv("name")))
        #print "parameter_set_file: %s" %(parameter_set_file)
        i_o.write_matrix(parameter_set, parameter_set_file)

        sy.setenv("PARAMETERSET", repr(parameter_set))

        # call the producer
        if self._objective_function == "Physical_Properties_Loss":
            producer_script = i_o.check_source_code(["parallel_jobs/parameter_substance_ensemble_temperature_variation.py"])[0]
        if self._objective_function == "QM_MM_Loss":
            producer_script = i_o.check_source_code(["parallel_jobs/parameter_variation.py"])[0]
        if self._objective_function == "PhysProp_QMMM_Loss":
            producer_script_Physical_Properties = i_o.check_source_code(["parallel_jobs/parameter_substance_ensemble_temperature_variation.py"])[0]
            producer_script_QM_MM = i_o.check_source_code(["parallel_jobs/parameter_variation.py"])[0]

        config_file = os.getenv("CONFIGFILE")
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            res = os.system("python %s %s %s" % (producer_script, config_file, parameter_set_file))
            if res != 0:
                print "ERROR in physical_properties_loss_new.get_function_values: Producer failed."
                tr.errorexit()
        if self._objective_function == "PhysProp_QMMM_Loss":
            res_QM_MM = subprocess.Popen(["python",os.path.join(self.__cwd,producer_script_QM_MM), config_file, parameter_set_file,"QM_MM"])
            ## wait for process to finish
            while 1:
                check_res = res_QM_MM.poll()
                time.sleep(1)
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.get_function_values: Producer for QM_MM failed."
                    tr.errorexit()

            res_Physical_Properties = subprocess.Popen(["python",os.path.join(self.__cwd,producer_script_Physical_Properties), config_file, parameter_set_file,"Physical_Properties"])
            ## wait for process to finish
            while 1:
                check_res = res_Physical_Properties.poll()
                time.sleep(1)
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.get_function_values: Producer for Physical_Properties failed."
                    tr.errorexit()



        if self._objective_function == "QM_MM_Loss" or self._objective_function == "PhysProp_QMMM_Loss":
            #print "len(parameter_set) = %s" %(repr(len(parameter_set)))
            #print "dim = %s" %(repr(self.get_dimension()))
            #self.define_working_directories_and_copy_config_file(len(parameter_set)-1, config_file)
            try:
                substance = self.__config.get("QM_MM","substance")
            except:
                print "ERROR in physical_properties_loss_new.get_function_values: You must indicate the name of your substance."
                tr.errorexit()

            if self._objective_function == "QM_MM_Loss":
                working_dir_file = os.path.join(os.getenv("GROW_HOME"), "working_directories.txt")
            if self._objective_function == "PhysProp_QMMM_Loss":
                working_dir_file = os.path.join(os.getenv("GROW_HOME_QM_MM"), "working_directories.txt")
            fwd = open(working_dir_file, "w")
          
            for k in range(len(parameter_set)):
                if self._objective_function == "QM_MM_Loss":
                    orig_working_dir = os.path.join(os.getenv("GROW_HOME"),os.getenv("name")+".%d" %(begin + k))
                if self._objective_function == "PhysProp_QMMM_Loss":
                    orig_working_dir = os.path.join(os.getenv("GROW_HOME_QM_MM"),os.getenv("name")+".%d" %(begin + k))
                working_dir = os.path.join(orig_working_dir,substance)
                try:
                    os.stat(working_dir)
                except:
                    os.makedirs(working_dir)

                os.system("cp %s %s" %(os.path.join(orig_working_dir,os.path.basename(config_file)),os.path.join(working_dir,os.path.basename(config_file))))
                fwd.write(working_dir+"\n")

            fwd.close()

        # call executer and collector
#        sy.setenv("PARAMETERSET", repr(parameter_set))

        # call the executer
        execution_scripts = i_o.check_source_code(["parallel_jobs/executer.py", "parallel_jobs/collector.py"])
        if self._objective_function == "Physical_Properties_Loss":
            simulation_script = i_o.check_source_code(["simulation/molecular_simulation.py"])[0]
        if self._objective_function == "QM_MM_Loss":
            simulation_script = i_o.check_source_code(["simulation/qm_mm.py"])[0]
        if self._objective_function == "PhysProp_QMMM_Loss":
            simulation_script_Physical_Properties = i_o.check_source_code(["simulation/molecular_simulation.py"])[0]
            simulation_script_QM_MM = i_o.check_source_code(["simulation/qm_mm.py"])[0]

        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            res = os.system("python %s %s %s" % (execution_scripts[0], config_file, simulation_script))
            if res != 0:
                print "ERROR in physical_properties_loss_new.get_function_values: Executer failed."
                tr.errorexit()
        if self._objective_function == "PhysProp_QMMM_Loss":
            res_QM_MM = subprocess.Popen(["python", os.path.join(self.__cwd,execution_scripts[0]), config_file, simulation_script_QM_MM, "QM_MM"])
            ## wait for subprocess to finish
            while 1:
                time.sleep(1)
                check_res = res_QM_MM.poll()
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.get_function_values: Executer for QM_MM failed."
                    tr.errorexit()

            res_Physical_Properties = subprocess.Popen(["python", os.path.join(self.__cwd,execution_scripts[0]), config_file, simulation_script_Physical_Properties, "Physical_Properties"])
            ## wait for subprocess to finish
            while 1:
                time.sleep(1)
                check_res = res_Physical_Properties.poll()
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.get_function_values: Executer for Physical_Properties failed."
                    tr.errorexit()

        
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "PhysProp_QMMM_Loss": 
            try:
                fits = self.__config.get("OPT","fits")
            except:
                print "ERROR in physical_properties_loss_new.get_function_values: You must indicate whether you want to perform properties fits or not."
                tr.errorexit()

            if fits != "y" and fits != "n":
                print "ERROR in physical_properties_loss_new.get_function_values: 'fits' is a boolean variable. Please choose between 'y' and 'n'. %s is not allowed." %(fits)
                tr.errorexit()
              
            # call the collector
            if self._objective_function == "Physical_Properties_Loss":
                res =  os.system("python %s %s %d %d %s" % (execution_scripts[1], os.path.join(os.getenv("GROW_HOME"), "property_dict.txt"), begin, begin + len(parameter_set) - 1, fits))
            if self._objective_function == "PhysProp_QMMM_Loss":
                res =  os.system("python %s %s %d %d %s %s" % (execution_scripts[1], os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), "property_dict.txt"), begin, begin + len(parameter_set) - 1, fits,"Physical_Properties"))
            if res != 0:
                print "ERROR in physical_properties_loss_new.get_function_values: Collector failed."
                tr.errorexit()
          
        loss_list = []

        if return_estimations:
            estimations_list = []

        if self._objective_function == "QM_MM_Loss" or self._objective_function == "Physical_Properties_Loss":
            for i in range(len(parameter_set)):
                properties_file = os.path.join(os.getenv("GROW_HOME"), os.getenv("name") + ".%d" %(begin + i), "properties.txt")
                #print "properties_file = %s" %(properties_file)
                self.set_estimations(properties_file)  
                loss_list.append(self.get_function_value(parameter_set[i]))
                if return_estimations:
                    #properties_list.append(self.get_estimations())
                    estimations_list.append(self.get_estimations())
        if self._objective_function == "PhysProp_QMMM_Loss":
            #for i in range(len(parameter_set)):
            #    properties_file = os.path.join(os.getenv("GROW_HOME_QM_MM"), os.getenv("name") + ".%d" %(begin + i), "properties.txt")
            #    self.set_estimations(properties_file)
            #    loss_list.append(self.get_function_value(parameter_set[i]))
            #    if return_estimations:
            #        #properties_list.append(self.get_estimations())
            #        estimations_list.append(self.get_estimations())
            #for i in range(len(parameter_set)):
            #    properties_file = os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), os.getenv("name") + ".%d" %(begin +i), "properties.txt")
            #    self.set_estimations(properties_file)
            #    loss_list.append(self.get_function_value(parameter_set[i]))
            #    if return_estimations:
            #        #properties_list.append(self.get_estimations())
            #        estimations_list.append(self.get_estimations())
            for i in range(len(parameter_set)):
                properties_file_Physical_Properties = os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), os.getenv("name") + ".%d" %(begin +i), "properties.txt")
                properties_file_QM_MM = os.path.join(os.getenv("GROW_HOME_QM_MM"), os.getenv("name") + ".%d" %(begin +i), "properties.txt")
                self.set_estimations(properties_file_QM_MM,QM=True)
                self.set_estimations(properties_file_Physical_Properties,PP=True)
                self._estimations=self._estimations_QM_MM
                for estimation_PP in self._estimations_Physical_Properties:
                    self._estimations.append(estimation_PP)
                loss_list.append(self.get_function_value(parameter_set[i]))
                if return_estimations:
                    #properties_list.append(self.get_estimations())
                    estimations_list.append(self.get_estimations())

        if return_estimations:
            #return [loss_list, properties_list]
            return [loss_list, estimations_list]
        else:
            print "loss_list:", loss_list 
            return loss_list


    def set_estimations(self, properties_file,QM=False,PP=False):
        """ sets estimations """
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            self._estimations = i_o.get_properties(properties_file)
        if self._objective_function == "PhysProp_QMMM_Loss":
            if QM == True:
                self._estimations_QM_MM = i_o.get_properties(properties_file)
            if PP == True:
                self._estimations_Physical_Properties = i_o.get_properties(properties_file)         

    def __change_parameter(self, y, h, i):
        """ help function: changes parameter x to x+h at index i """
        y[i] = y[i] + h

        return y


    def set_estimation_gradients(self, x):
        """ sets gradients of estimations using finite differences """
        try:
            h = self.__config.get("OPT", "h")
        except:
            print "ERROR in physical_properties_loss_new.set_estimation_gradients: You must indicate the step length h for the approximation of the derivatives."
            tr.errorexit()

        h = string.split(h)
        h = su.float_vector(h, check=True)
        if type(h) != list:
            print "ERROR in physical_properties_loss_new.set_estimation_gradients: Your indicated value for 'h' ('%s') is not numeric." % (h)
            tr.errorexit()

        dim = self.get_dimension()
        if len(h) < dim:
            for k in range(dim - 1):
                h.append(h[0])

        if int(os.getenv("loop")) == 1:
            # original parameter has to be simulated
            parameter_set = [x]
            begin = 0
        else:
            # original parameter was already simulated
            parameter_set = []
            begin = 1

        for k in range(dim):
            y = list(x)
            p = self.__change_parameter(y, h[k], k)
            parameter_set.append(p)

        sy.setenv("BEGIN", repr(begin))
         
        ## slot for the computation of physical properties ##
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            parameter_set_file = os.path.join(os.getenv("GROW_HOME"), "parameter_set.%s" % (os.getenv("name")))
            i_o.write_matrix(parameter_set, parameter_set_file)
        if self._objective_function == "PhysProp_QMMM_Loss":
            # QM_MM
            #print "in physical_properties_loss_new.py: os.getenv(GROW_HOME_QM_MM):",os.getenv("GROW_HOME_QM_MM")
            parameter_set_file_QM_MM = os.path.join(os.getenv("GROW_HOME_QM_MM"), "parameter_set.%s" % (os.getenv("name")))
            i_o.write_matrix(parameter_set, parameter_set_file_QM_MM)
            # Physical_Properties
            parameter_set_file_Physical_Properties = os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), "parameter_set.%s" % (os.getenv("name")))
            i_o.write_matrix(parameter_set, parameter_set_file_Physical_Properties)

        # call the producer
        if self._objective_function == "Physical_Properties_Loss":
            producer_script = i_o.check_source_code(["parallel_jobs/parameter_substance_ensemble_temperature_variation.py"])[0]
        if self._objective_function == "QM_MM_Loss":
            producer_script = i_o.check_source_code(["parallel_jobs/parameter_variation.py"])[0]
        if self._objective_function == "PhysProp_QMMM_Loss":
            producer_script_Physical_Properties = i_o.check_source_code(["parallel_jobs/parameter_substance_ensemble_temperature_variation.py"])[0]
            producer_script_QM_MM = i_o.check_source_code(["parallel_jobs/parameter_variation.py"])[0]

        config_file = os.getenv("CONFIGFILE")
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            res = os.system("python %s %s %s" % (producer_script, config_file, parameter_set_file))
            if res != 0:
                print "ERROR in physical_properties_loss_new.set_estimation_gradients: Producer failed."
                tr.errorexit()
        if self._objective_function == "PhysProp_QMMM_Loss":
            #print "pyhsical_properties_loss_new.py: [\"python\",os.path.join(self.__cwd,producer_script_QM_MM), config_file, parameter_set_file_QM_MM]:",["python",os.path.join(self.__cwd,producer_script_QM_MM), config_file, parameter_set_file_QM_MM]
            res_QM_MM = subprocess.Popen(["python",os.path.join(self.__cwd,producer_script_QM_MM), config_file, parameter_set_file_QM_MM, "QM_MM"])
            ## wait for process to finish
            while 1:
                check_res = res_QM_MM.poll()
                time.sleep(1)
                #print "pyhsical_properties_loss_new.py: check_res:",check_res
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.set_estimation_gradient: Producer for QM_MM failed."
                    tr.errorexit()
            #print "pyhsical_properties_loss_new.py: [\"python\",os.path.join(self.__cwd,producer_script_Physical_Properties), config_file, parameter_set_file_Physical_Properties]:",["python",os.path.join(self.__cwd,producer_script_Physical_Properties), config_file, parameter_set_file_Physical_Properties, "Physical_Properties"]
            res_Physical_Properties = subprocess.Popen(["python",os.path.join(self.__cwd,producer_script_Physical_Properties), config_file, parameter_set_file_Physical_Properties, "Physical_Properties"])
            ## wait for process to finish
            while 1:
                check_res = res_Physical_Properties.poll()
                time.sleep(1)
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.set_estimation_gradient: Producer for Physical_Properties failed."
                    tr.errorexit()

        # defining the working directories and copy config files into them
        if self._objective_function == "QM_MM_Loss" or self._objective_function == "PhysProp_QMMM_Loss":
            #print "dim = %s" %(repr(dim))
            #self.define_working_directories_and_copy_config_file(dim, config_file)
            try:
                substance = self.__config.get("QM_MM","substance")
            except:
                try:
                    substance = self.__config.get("OPT","substance")
                except:
                    print "ERROR in physical_properties_loss_new.set_estimation_gradients: You must indicate the name of your substance in [QM_MM] or [OPT]."
                    tr.errorexit()

            if self._objective_function == "QM_MM_Loss":
                working_dir_file = os.path.join(os.getenv("GROW_HOME"), "working_directories.txt")
            if self._objective_function == "PhysProp_QMMM_Loss":
                working_dir_file = os.path.join(os.getenv("GROW_HOME_QM_MM"), "working_directories.txt")

            fwd = open(working_dir_file, "w")
            for k in range(dim + 1):
                if self._objective_function == "QM_MM_Loss":
                    orig_working_dir = os.path.join(os.getenv("GROW_HOME"),os.getenv("name")+".%d" %(k))
                if self._objective_function == "PhysProp_QMMM_Loss":
                    orig_working_dir = os.path.join(os.getenv("GROW_HOME_QM_MM"),os.getenv("name")+".%d" %(k))

                working_dir = os.path.join(orig_working_dir,substance)
                try:
                    os.stat(working_dir)
                except:
                    os.makedirs(working_dir)
                os.system("cp %s %s" %(os.path.join(orig_working_dir,os.path.basename(config_file)),os.path.join(working_dir,os.path.basename(config_file))))
                fwd.write(working_dir+"\n")

            fwd.close()


        # call executer and collector
        sy.setenv("PARAMETERSET", repr(parameter_set))
        execution_scripts = i_o.check_source_code(["parallel_jobs/executer.py", "parallel_jobs/collector.py"])
        if self._objective_function == "Physical_Properties_Loss":
            simulation_script = i_o.check_source_code(["simulation/molecular_simulation.py"])[0]
        if self._objective_function == "QM_MM_Loss":
            simulation_script = i_o.check_source_code(["simulation/qm_mm.py"])[0]
        if self._objective_function == "PhysProp_QMMM_Loss":
            simulation_script_Physical_Properties = i_o.check_source_code(["simulation/molecular_simulation.py"])[0]
            simulation_script_QM_MM = i_o.check_source_code(["simulation/qm_mm.py"])[0]

        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            res = os.system("python %s %s %s" % (execution_scripts[0], config_file, simulation_script))
            if res != 0:
                print "ERROR in physical_properties_loss_new.pset_estimation_gradients: Executer failed."
                print "execution_scripts[0] = ",execution_scripts[0]
                print "config_file = ",config_file
                print "simulation_script = ",simulation_script
                tr.errorexit()
        if self._objective_function == "PhysProp_QMMM_Loss":
            res_QM_MM = subprocess.Popen(["python", os.path.join(self.__cwd,execution_scripts[0]), config_file, simulation_script_QM_MM, "QM_MM"])
            ## wait for process to finish
            while 1:
                check_res = res_QM_MM.poll()
                time.sleep(1)
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.set_estimation_gradient: Executer for QM_MM failed."
                    print "execution_scripts[0] = ",execution_scripts[0]
                    print "config_file = ",config_file
                    print "simulation_script = ",simulation_script_QM_MM
                    tr.errorexit()

            res_Physical_Properties = subprocess.Popen(["python", os.path.join(self.__cwd,execution_scripts[0]), config_file, simulation_script_Physical_Properties, "Physical_Properties"])
            ## wait for process to finish
            while 1:
                check_res = res_Physical_Properties.poll()
                time.sleep(1)
                if check_res == 0:
                    break
                if check_res == 1:
                    print "ERROR in physical_properties_loss_new.set_estimation_gradient: Executer for Physical_Properties failed."
                    print "execution_scripts[0] = ",execution_scripts[0]
                    print "config_file = ",config_file
                    print "simulation_script = ",simulation_script_Physical_Properties
                    tr.errorexit()


        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "PhysProp_QMMM_Loss":          
            try:
                fits = self.__config.get("OPT","fits")
            except:
                print "ERROR in physical_properties_loss_new.set_estimation_gradients: You must indicate whether you want to perform properties fits or not."
                tr.errorexit()

            if fits != "y" and fits != "n":
                print "ERROR in physical_properties_loss_new.set_estimation_gradients: 'fits' is a boolean variable. Please choose between 'y' and 'n'. %s is not allowed." %(fits)
                tr.errorexit()
              
            # collector
            if self._objective_function == "Physical_Properties_Loss":
                res =  os.system("python %s %s %d %d %s" % (execution_scripts[1], os.path.join(os.getenv("GROW_HOME"), "property_dict.txt"), begin, dim, fits))
            if self._objective_function == "PhysProp_QMMM_Loss":
                res =  os.system("python %s %s %d %d %s %s" % (execution_scripts[1], os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), "property_dict.txt"), begin, dim, fits, "Physical_Properties"))
            if res != 0:
                print "ERROR in physical_properties_loss_new.set_estimation_gradients: Collector failed."
                tr.errorexit()

        # original x
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            properties_file0 = os.path.join(os.getenv("GROW_HOME"), os.getenv("name") + ".0", "properties.txt")
            self.set_estimations(properties_file0)
            properties_0 = list(self._estimations)
            ## end slot ##
        if self._objective_function == "PhysProp_QMMM_Loss":
            properties_file0_QM_MM = os.path.join(os.getenv("GROW_HOME_QM_MM"), os.getenv("name") + ".0", "properties.txt")
            self.set_estimations(properties_file0_QM_MM,QM=True)
            properties_0_QM_MM = list(self._estimations_QM_MM)

            properties_file0_Physical_Properties = os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), os.getenv("name") + ".0", "properties.txt")
            self.set_estimations(properties_file0_Physical_Properties,PP=True)
            properties_0_Physical_Properties = list(self._estimations_Physical_Properties)

            properties_0 = properties_0_QM_MM
            for prop0_PP in properties_0_Physical_Properties:
                properties_0.append(prop0_PP)
          
        # initialize gradients
        #if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
        self._gradients_of_estimations = []
        for i in range(len(properties_0)):
            self._gradients_of_estimations.append([])
            for k in range(dim):
                self._gradients_of_estimations[i].append(0.0)
        # probably not needed
        #if self._objective_function == "PhysProp_QMMM_Loss":
        #    self._gradients_of_estimations_QM_MM = []
        #    self._gradients_of_estimations_Physical_Properties = []
        #    for i in range(len(properties_0_QM_MM)):
        #        self._gradients_of_estimations_QM_MM.append([])
        #        for k in range(dim):
        #            self._gradients_of_estimations_QM_MM[i].append(0.0)
        #    for i in range(len(properties_0_Physical_Properties)):
        #        self._gradients_of_estimations_Physical_Properties.append([])
        #        for k in range(dim):
        #            self._gradients_of_estimations_Physical_Properties[i].append(0.0)

        #properties_tmp=[]  
        for k in range(1, dim + 1):
            if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
                properties_file = os.path.join(os.getenv("GROW_HOME"), os.getenv("name") + ".%d" % (k), "properties.txt")
                #print "properties_file = %s" %(properties_file)
                self.set_estimations(properties_file)
                properties_k = list(self._estimations)
                ### check if props can be appended via .append ### seems to work - not tested yet
                #print "properties_k:",properties_k
                #for prop in list(self._estimations):
                #    properties_k.append(prop)
                #print "properties_k_appended:",properties_k

                if len(properties_k) != len(properties_0):
                    print "ERROR in physical_properties_loss_new.set_estimation_gradients: The length of the %dth property vector (%d) does not correspond to the length of the original one (%d)" % (k, len(properties_k), len(properties_0))
                    tr.errorexit()

            if self._objective_function == "PhysProp_QMMM_Loss":
                properties_file_QM_MM = os.path.join(os.getenv("GROW_HOME_QM_MM"), os.getenv("name") + ".%d" % (k), "properties.txt")
                self.set_estimations(properties_file_QM_MM,QM=True)
                properties_k_QM_MM = list(self._estimations_QM_MM)

                properties_file_Physical_Properties = os.path.join(os.getenv("GROW_HOME_PHYSICAL_PROPERTIES"), os.getenv("name") + ".%d" % (k), "properties.txt")
                self.set_estimations(properties_file_Physical_Properties,PP=True)
                properties_k_Physical_Properties = list(self._estimations_Physical_Properties)

                properties_k = properties_k_QM_MM
                for prop_PP in properties_k_Physical_Properties:
                    properties_k.append(prop_PP)
                #properties_tmp.append(properties_k)

                if len(properties_k_QM_MM) != len(properties_0_QM_MM):
                    print "ERROR in physical_properties_loss_new.set_estimation_gradients: The length of the %dth property vector (QM) (%d) does not correspond to the length of the original one (%d)" % (k, len(properties_k_QM_MM), len(properties_0_QM_MM))
                    tr.errorexit()

                if len(properties_k_Physical_Properties) != len(properties_0_Physical_Properties):
                    print "ERROR in physical_properties_loss_new.set_estimation_gradients: The length of the %dth property vector (PP) (%d) does not correspond to the length of the original one (%d)" % (k, len(properties_k_Physical_Properties), len(properties_0_Physical_Properties))
                    tr.errorexit()
        
            finite_differences = m.mult_vector(1. / h[k-1], m.vector_sum(properties_k, m.mult_vector(-1., properties_0)))
            for i in range(len(properties_0)):
                self._gradients_of_estimations[i][k-1] = finite_differences[i]

        #test_file = open('/home/rstric2s/properties.txt','w')
        #test_file.write("properties_0 \t properties_1 \t properties_2 \t properties_3 \t properties_4 \n")
        #for i in range(len(properties_0)):
        #    tmp_props = str(properties_0[i]) + "\t"
        #    for j in range(len(properties_tmp)):
        #        tmp_props = tmp_props + str(properties_tmp[j][i]) + '\t'
        #    test_file.write(tmp_props + "\n")
        #test_file.close()

        #reset properties
        if self._objective_function == "Physical_Properties_Loss" or self._objective_function == "QM_MM_Loss":
            self.set_estimations(properties_file0)
        if self._objective_function == "PhysProp_QMMM_Loss":
            self.set_estimations(properties_file0_QM_MM,QM=True)
            self.set_estimations(properties_file0_Physical_Properties,PP=True)
            self._estimations = self._estimations_QM_MM
            for estimation_PP in self._estimations_Physical_Properties:
                self._estimations.append(estimation_PP)

    def set_estimation_hessians(self, x):
        """ sets hessians of estimations using finite differences """
        pass

        
    def get_function_value(self, x):
        """ return the loss function value of a vector x """
        return Loss.get_function_value(self, x)


    def get_gradient(self, x):
        """ return the gradient of the loss function at x """
        #print "in physical_properties_loss_new.get_gradient: x=",x
        self.set_estimation_gradients(x)
        gradient = Loss.get_gradient(self, x)

        # write summary files
        if self._objective_function == "Physical_Properties_Loss":
            opt.complete_summary_file(int(os.getenv("loop")) - 1, x, gradient, m.norm(gradient), self.get_estimations(), self.get_targets(), self.get_function_value(x))
            opt.complete_summary_table(int(os.getenv("loop")) - 1, self.get_estimations(), self.get_targets(), self.get_function_value(x))
        if self._objective_function == "QM_MM_Loss":
            i_o.complete_summary_file(int(os.getenv("loop")) - 1, x, gradient, m.norm(gradient), [], self.get_targets(), self.get_function_value(x))
            i_o.complete_summary_table(int(os.getenv("loop")) - 1, [], self.get_targets(), self.get_function_value(x))
        if self._objective_function == "PhysProp_QMMM_Loss":
            # Physical Properties
            #opt.complete_summary_file(int(os.getenv("loop")) - 1, x, gradient, m.norm(gradient), self.get_estimations(), self.get_targets(), self.get_function_value(x))
            #opt.complete_summary_table(int(os.getenv("loop")) - 1, self.get_estimations(), self.get_targets(), self.get_function_value(x))
            # QM MM
            i_o.complete_summary_file(int(os.getenv("loop")) - 1, x, gradient, m.norm(gradient), [], self.get_targets(), self.get_function_value(x),hybrid_opt=True)
            i_o.complete_summary_table(int(os.getenv("loop")) - 1, [], self.get_targets(), self.get_function_value(x)),
          
        return gradient


    def get_hessian(self, x):
        """ return the Hessian of the loss function at x """
        self.set_estimation_hessians(x)
          
        return Loss.get_hessian(self, x)

# unused function due to determination of Physical_Properties_Loss or QM_MM_Loss in the same script.
#    def define_working_directories_and_copy_config_file(self, dim, config_file):
#        """ defining the working directories and copy config files into them """
#        try:
#            substance = self.__config.get("QM_MM","substance")
#        except:
#            print "ERROR in physical_properties_loss_new.get_hessian: You must indicate the name of your substance."
#            tr.errorexit()
#
#        working_dir_file = os.path.join(os.getenv("GROW_HOME"), "working_directories.txt")
#        fwd = open(working_dir_file, "w")
#        for k in range(dim + 1):
#            orig_working_dir = os.path.join(os.getenv("GROW_HOME"),os.getenv("name")+".%d" %(k))
#            working_dir = os.path.join(orig_working_dir,substance)
#            try:
#                os.stat(working_dir)
#            except:
#                os.makedirs(working_dir)
#
#            os.system("cp %s %s" %(os.path.join(orig_working_dir,os.path.basename(config_file)),os.path.join(working_dir,os.path.basename(config_file))))
#            fwd.write(working_dir+"\n")
#
#        fwd.close()


    def __del__(self):
        """ destructor """
        del self
