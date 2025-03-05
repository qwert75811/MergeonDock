# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:55:07 2024

@author: Xhamrock Studio
"""

import os
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QMessageBox
from PyQt5.QtCore import Qt

import subprocess
import shutil



class Dock_setting(QtCore.QObject):
    def __init__(self, ui, all_parameters, analysis_basic_tab):
        super().__init__()
        self.ui = ui
        self.all_parameters = all_parameters
        self.analysis_basic_tab = analysis_basic_tab
        self.ui.pushButton_dockingbutton.clicked.connect(self.run_docking)
        self.dock_process = None
        self.current_index = 0
        self.total_ligands = 0
        self.ui.progressBar.setValue(0)
        self.is_docking = False
        
    def run_docking(self):
        self.current_index = 0
        self.set_progress_value(0)
        self.output_log_path = []
        self.output_file_path = []
        
        self.total_ligands = len(self.all_parameters.output_prepared_ligands_path)
        if self.total_ligands > 0:
            self.is_docking = True
            self.ui.pushButton_dockingbutton.setText("Stop")
            self.ui.pushButton_dockingbutton.clicked.disconnect()
            self.ui.pushButton_dockingbutton.clicked.connect(self.cancel_process)
            self.run_next_process()
        else:
            error_message = "Warning message:"
            error_window = QtWidgets.QMessageBox()
            error_window.setIcon(QtWidgets.QMessageBox.Critical)
            error_window.setWindowTitle("No ligands found")
            error_window.setText("No ligands are input, please upload ligands again.")
            error_window.setInformativeText(error_message)
            error_window.setStandardButtons(QtWidgets.QMessageBox.Ok)
            error_window.exec_()  
       
    def run_next_process(self):
        if self.current_index < self.total_ligands and self.is_docking:
            ligand = self.all_parameters.output_prepared_ligands_path[self.current_index]
            ligand_filename = os.path.basename(ligand)
            ligand_name = os.path.splitext(ligand_filename)[0]
           
            output_folder_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{ligand_name}"))
            if os.path.exists(output_folder_path):
                shutil.rmtree(output_folder_path)
                os.mkdir(output_folder_path)
            else:
                os.mkdir(output_folder_path)
                
            output_path = os.path.normpath(os.path.join(output_folder_path, f"{self.all_parameters.output_prepared_receptor_name}_{ligand_name}.pdbqt"))
            self.output_log = os.path.normpath(os.path.join(output_folder_path, f"{self.all_parameters.output_prepared_receptor_name}_{ligand_name}_log.txt"))
            self.output_log_path.append(self.output_log)
            self.output_file_path.append(output_path)
            
            
            if self.all_parameters.scoring_function == "vina":
                vina_command = (f'"{self.all_parameters.vina_path}" '
                                f'--receptor "{self.all_parameters.output_prepared_receptor_path}" '
                                f'--ligand "{ligand}" '
                                f'--scoring vina '
                                f'--center_x {self.all_parameters.gridcenter_X} '
                                f'--center_y {self.all_parameters.gridcenter_Y} '
                                f'--center_z {self.all_parameters.gridcenter_Z} '
                                f'--size_x {self.all_parameters.gridsize_X} '
                                f'--size_y {self.all_parameters.gridsize_Y} '
                                f'--size_z {self.all_parameters.gridsize_Z} '
                                f'--cpu {self.all_parameters.cpu_value} '
                                f'--seed {self.all_parameters.seed_value} '
                                f'--exhaustiveness {self.all_parameters.exhaustiveness} '
                                f'--num_modes {self.all_parameters.poses} '
                                f'--spacing {self.all_parameters.spacing_value} '
                                f'--verbosity {self.all_parameters.verbosity_value} '
                                f'--out "{output_path}"'
                )
                self.start_process(vina_command, self.current_index + 1, self.total_ligands)
                
                
                
            elif self.all_parameters.scoring_function == "autodock4":
                work_receptor_path = os.path.normpath(os.path.join(output_folder_path, f"{self.all_parameters.output_prepared_receptor_name}.pdbqt"))
                work_ligand_path = os.path.normpath(os.path.join(output_folder_path, f"{ligand_filename}"))
                
                shutil.copy(f"{self.all_parameters.output_prepared_receptor_path}", work_receptor_path)
                shutil.copy(f"{ligand}", work_ligand_path)
 
                out_gpf_path = os.path.normpath(os.path.join(output_folder_path, f"{self.all_parameters.output_prepared_receptor_name}_{ligand_name}.gpf"))
                out_glg_path = os.path.normpath(os.path.join(output_folder_path, f"{self.all_parameters.output_prepared_receptor_name}_{ligand_name}.glg"))
                
                gpf_command = (f'{self.all_parameters.autodock4_run_prepare_gpf} '
                               f'-l "{ligand}" -r "{self.all_parameters.output_prepared_receptor_path}" '
                               f'-o "{out_gpf_path}" '
                               f'-p npts={self.all_parameters.gridsize_X},{self.all_parameters.gridsize_Y},{self.all_parameters.gridsize_Z} '
                               f'-p gridcenter={self.all_parameters.gridcenter_X},{self.all_parameters.gridcenter_Y},{self.all_parameters.gridcenter_Z}'
                )
                subprocess.run(gpf_command)
                
                
                autogrid4_command = f'"{self.all_parameters.autodock4_autogrid4_path}" -p "{out_gpf_path}" -l "{out_glg_path}"'
                
                
                vina_command = (f'"{self.all_parameters.vina_path}" ' 
                                f'--ligand "{work_ligand_path}" '
                                f'--maps {self.all_parameters.output_prepared_receptor_name} '
                                f'--scoring ad4 '
                                f'--cpu {self.all_parameters.cpu_value} '
                                f'--seed {self.all_parameters.seed_value} '
                                f'--exhaustiveness {self.all_parameters.exhaustiveness} '
                                f'--num_modes {self.all_parameters.poses} '
                                f'--spacing {self.all_parameters.spacing_value} '
                                f'--verbosity {self.all_parameters.verbosity_value} '
                                f'--out "{output_path}"'
                )
                
                os.chdir(output_folder_path)
                self.grid_process = QtCore.QProcess()
                self.grid_process.finished.connect(lambda: self.start_process(vina_command, self.current_index + 1, self.total_ligands))
                self.grid_process.start(autogrid4_command)
                
    
    
    

    def start_process(self, command, current, total):
        self.dock_process = QtCore.QProcess()
        self.dock_process.finished.connect(lambda: self.dock_process_finished(current, total))
        self.dock_process.start(command)
        
        

    def dock_process_finished(self, current, total):
        stdout_data = self.dock_process.readAllStandardOutput().data().decode()
        stderr_data = self.dock_process.readAllStandardError().data().decode()
        with open(self.output_log, "w") as file:
            file.write(stdout_data)
        
        folder_path = os.path.dirname(self.output_log)
        error_path = os.path.normpath(os.path.join(folder_path, "error_log.txt"))
        if stderr_data.strip():
            with open(error_path, "w") as file:
                file.write(stderr_data)
        print(stdout_data)
        print(stderr_data)

        progress_value = int((current / total) * 100)
        self.set_progress_value(progress_value)
        
        if current == total:
            self.set_progress_value(100)
            self.ui.pushButton_dockingbutton.setText("Docking")
            self.ui.pushButton_dockingbutton.clicked.disconnect()
            self.ui.pushButton_dockingbutton.clicked.connect(self.run_docking)
            self.is_docking = False
            os.chdir(self.all_parameters.work_directory)
            self.work_log()
        else:
            self.current_index += 1
            self.run_next_process()
        
        

    def set_progress_value(self, value):
        self.ui.progressBar.setValue(value)

    
    def cancel_process(self):
        if self.dock_process:
            self.dock_process.kill()
        self.ui.pushButton_dockingbutton.setText("Docking")
        self.set_progress_value(0)
        self.current_index = 0
        self.ui.pushButton_dockingbutton.clicked.disconnect()
        self.ui.pushButton_dockingbutton.clicked.connect(self.run_docking)
        self.is_docking = False
        os.chdir(self.all_parameters.work_directory)
        
        if self.dock_process:
            self.dock_process.disconnect()
            self.dock_process = None
    
    
    
    def work_log(self):
        output_filename_list = []
        for output_file_name in self.output_file_path:
            output_filename = os.path.basename(output_file_name)
            output_filename_list.append(output_filename)
        
        output_logname_list = []
        for output_log_name in self.output_log_path:
            output_logname = os.path.basename(output_log_name)
            output_logname_list.append(output_logname)
            
        
        
        reflig_basename = os.path.splitext(os.path.basename(self.all_parameters.ref_ligand_picked_path))[0]
        
        
        work_log_text = f"""
Created by MergeonDock
        
Work directory: {self.all_parameters.work_directory}
        
Receptor: {self.all_parameters.output_prepared_receptor_name}
Ligands: {self.all_parameters.output_prepared_ligands_name}
Ref ligand: {reflig_basename}
        
Scoring function: {self.all_parameters.scoring_function}
        
output files: {output_filename_list}

output logs: {output_logname_list}
        
"""
        log_path = os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.output_prepared_receptor_name}.cdl")
        
        with open(log_path, "w") as file:
            file.write(work_log_text)
            
        self.all_parameters.cdl_path = log_path
        self.ask_for_analysis()
    
    
    
    def ask_for_analysis(self):
        question_window = QMessageBox()
        question_window.setIcon(QMessageBox.Question)
        question_window.setWindowTitle("Completed")
        question_window.setText("Docking completed. Do you want to move to Analysis section?")
        question_window.setStandardButtons(QMessageBox.Yes | QMessageBox.No)  # 添加 Yes 和 No 按鈕
        question_window.setDefaultButton(QMessageBox.No)  # 預設選擇 No 按鈕
        
        question_window.setWindowModality(Qt.ApplicationModal)
        question_window_reply = question_window.exec_()
        
        if question_window_reply == question_window.Yes:
            self.analysis_basic_tab.auto_load_from_dock_tab(self.all_parameters.cdl_path)
        else:
            print("cancel")
        

                
        
