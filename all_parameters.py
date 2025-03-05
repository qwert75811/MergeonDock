# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:42:01 2024

@author: Xhamrock Studio
"""

import os
import sys

class Parameters_storage:
    def __init__(self):
        if getattr(sys, 'frozen', False):
            self.install_path = os.path.dirname(sys.executable)
                                                       
            
            
            self.autodock4_folder_path = os.path.normpath(os.path.join(self.install_path,"_internal", "autodock4"))
            self.python2_path = os.path.normpath(os.path.join(self.autodock4_folder_path, "python.exe"))
            self.autodock4_script_path = os.path.normpath(os.path.join(self.autodock4_folder_path, "Scripts"))
            self.autodock4_autogrid4_path = os.path.normpath(os.path.join(self.autodock4_script_path, "autogrid4.exe"))
            
            self.vina_path = os.path.normpath(os.path.join(self.autodock4_script_path, "vina_1.2.5_win.exe"))
  
            self.autodock4_Utilities24_path = os.path.normpath(os.path.join(self.autodock4_folder_path, "Lib", "site-packages", "AutoDockTools", "Utilities24"))
            
            self.autodock4_prepare_receptor_path = os.path.normpath(os.path.join(self.autodock4_Utilities24_path, "prepare_receptor4.py"))
            self.autodock4_prepare_ligand4_path = os.path.normpath(os.path.join(self.autodock4_Utilities24_path, "prepare_ligand4.py"))
            self.autodock4_prepare_gpf_path = os.path.normpath(os.path.join(self.autodock4_Utilities24_path, "prepare_gpf.py"))
            
            self.autodock4_run_prepare_receptor = f'"{self.python2_path}" "{self.autodock4_prepare_receptor_path}"' 
            self.autodock4_run_prepare_ligands = f'"{self.python2_path}" "{self.autodock4_prepare_ligand4_path}"' 
            self.autodock4_run_prepare_gpf = f'"{self.python2_path}" "{self.autodock4_prepare_gpf_path}"' 
            
        else:
            self.install_path = os.path.dirname(os.path.abspath(__file__))
            self.site_package_path = os.path.dirname(self.install_path)
            self.lib_path = os.path.dirname(self.site_package_path)
            self.env_path = os.path.dirname(self.lib_path)
            self.script_path = os.path.normpath(os.path.join(self.env_path, "Scripts"))
            
            self.autodock4_script_path = os.path.normpath(os.path.join(self.site_package_path, "autodock4", "Scripts"))
            self.autodock4_autogrid4_path = os.path.normpath(os.path.join(self.autodock4_script_path, "autogrid4.exe"))
            
            self.vina_path = os.path.normpath(os.path.join(self.autodock4_script_path, "vina_1.2.5_win.exe"))
            
            
            
            
            self.python2_path = os.path.normpath(os.path.join(self.site_package_path, "autodock4", "python.exe"))
            self.autodock4_Utilities24_path = os.path.normpath(os.path.join(self.site_package_path, "autodock4", "Lib", "site-packages", "AutoDockTools", "Utilities24"))
            
            self.autodock4_prepare_receptor_path = os.path.normpath(os.path.join(self.autodock4_Utilities24_path, "prepare_receptor4.py"))
            self.autodock4_prepare_ligand4_path = os.path.normpath(os.path.join(self.autodock4_Utilities24_path, "prepare_ligand4.py"))
            self.autodock4_prepare_gpf_path = os.path.normpath(os.path.join(self.autodock4_Utilities24_path, "prepare_gpf.py"))
            
            self.autodock4_run_prepare_receptor = f"{self.python2_path} {self.autodock4_prepare_receptor_path}" 
            self.autodock4_run_prepare_ligands = f"{self.python2_path} {self.autodock4_prepare_ligand4_path}" 
            self.autodock4_run_prepare_gpf = f"{self.python2_path} {self.autodock4_prepare_gpf_path}" 
        

    

        
        self.autodock4_prepare_gpf_auto_center_command = "-y" 
        
        self.autodock_prepare_receptor_custom_command = ""
        self.autodock_prepare_ligands_custom_command = ""
        
        
        self.work_directory = None
        
          
        self.receptor_prepare_method = "ad4"
        self.receptor_prepare_opt_switch = False
        
        self.receptor_opt_parameters_keys = ["rec_ad4_A", 
                                             "rec_ad4_A_combobox", 
                                             "rec_ad4_C", 
                                             "rec_ad4_e", 
                                             "rec_ad4_w", 
                                             "rec_ad4_p", 
                                             "rec_ad4_p_lineedit", 
                                             "rec_ad4_d", 
                                             "rec_ad4_d_lineedit", 
                                             "rec_ad4_U",
                                             "rec_ad4_U_nphs",
                                             "rec_ad4_U_lps",
                                             "rec_ad4_U_waters",
                                             "rec_ad4_U_nonstdres"
                                             ]
        self.receptor_opt_parameters_dict = {key: "" for key in self.receptor_opt_parameters_keys}
        
        
        
        self.ligands_prepare_method = "ad4"
        self.ligands_prepare_opt_switch = False
        
        self.ligands_opt_parameters_keys = ["lig_ad4_A", 
                                             "lig_ad4_A_combobox", 
                                             "lig_ad4_C", 
                                             "lig_ad4_Z", 
                                             "lig_ad4_g",
                                             "lig_ad4_s",
                                             "lig_ad4_w",
                                             "lig_ad4_F",
                                             "lig_ad4_R",
                                             "lig_ad4_R_lineedit",
                                             "lig_ad4_I",
                                             "lig_ad4_I_lineedit",
                                             "lig_ad4_p", 
                                             "lig_ad4_p_lineedit", 
                                             "lig_ad4_d", 
                                             "lig_ad4_d_lineedit", 
                                             "lig_ad4_U",
                                             "lig_ad4_U_nphs",
                                             "lig_ad4_U_lps",
                                             "lig_ad4_B",
                                             "lig_ad4_B_backbone",
                                             "lig_ad4_B_amide",
                                             "lig_ad4_B_guanidinium",
                                             ]
        self.ligands_opt_parameters_dict = {key: "" for key in self.ligands_opt_parameters_keys}
        
        
 
        
        
        
        self.scoring_function = "vina"
        
        
        self.input_ligands_path = []
        self.input_ligands_name = []
        
        self.input_receptor_path = None
        self.input_receptor_name = None
        
        self.input_flexible = None
         
        self.output_prepared_receptor_path = ""
        self.output_prepared_receptor_name = ""
        
        self.output_prepared_ligands_path = []
        self.output_prepared_ligands_name = []
        
        self.ref_prepared_ligands_path = []
        self.ref_prepared_ligands_name = []
        self.ref_ligand_picked_path = ""
          
        self.gridcenter_X = 0
        self.gridcenter_Y = 0
        self.gridcenter_Z = 0
        self.gridsize_X = 40
        self.gridsize_Y = 40
        self.gridsize_Z = 40
        self.exhaustiveness = 32
        self.poses = 2
        self.seed_value = 0
        self.spacing_value = 0.375
        self.verbosity_value = 2
        self.cpu_value = 0
        
        self.cdl_path = ""
        
