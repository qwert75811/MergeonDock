# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 11:09:41 2024

@author: Xhamrock Studio
"""

import os
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox


class Directory_setup():
    def __init__(self, ui, pymol_process, all_parameters):
        self.ui = ui
        self.all_parameters = all_parameters
        self.pymol_process = pymol_process
        
        
        self.ui.pushButton_setupworkdirectory.clicked.connect(self.button_setup_workdirectory)
        
    
    def button_setup_workdirectory(self):
        # 如果已經有設定工作目錄，詢問使用者是否要清空並重新選擇
        if self.all_parameters.work_directory:
            reply = QMessageBox.question(
                None,
                "Reset Work Directory",
                "A work directory is already set.\nDo you want to start a new one?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )

            if reply == QMessageBox.No:
                return  # 如果選擇 "No"，則不變更工作目錄

            # 選擇 "Yes"，則清空所有相關變數
            
        
        
        work_folder = QtWidgets.QFileDialog.getExistingDirectory()

        if work_folder:
            self.all_parameters.work_directory = work_folder
            self.clear_all_parameters()
            self.show_work_directory()
            self.button_activate()
        
        
            
    def show_work_directory(self):
        self.ui.textBrowser_workdirectory.clear()
        if self.all_parameters.work_directory:
            self.ui.textBrowser_workdirectory.append(self.all_parameters.work_directory)
        
        
    
    def button_activate(self):
        self.ui.pushButton_uploadreceptor.setEnabled(True)
        self.ui.pushButton_uploadligands.setEnabled(True)
        #self.ui.pushButton_uploadflexible.setEnabled(True)
        #self.ui.pushButton_uploadflexible.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(0, 41, 125)")
    
    
    def clear_all_parameters(self):
        """清空所有參數，模擬初始化狀態"""
        self.all_parameters.input_ligands_path = []
        self.all_parameters.input_ligands_name = []
        
        self.all_parameters.input_receptor_path = None
        self.all_parameters.input_receptor_name = None
        
        self.all_parameters.input_flexible = None
         
        self.all_parameters.output_prepared_receptor_path = ""
        self.all_parameters.output_prepared_receptor_name = ""
        
        self.all_parameters.output_prepared_ligands_path = []
        self.all_parameters.output_prepared_ligands_name = []
        
        self.all_parameters.ref_prepared_ligands_path = []
        self.all_parameters.ref_prepared_ligands_name = []
        self.all_parameters.ref_ligand_picked_path = ""
        
        self.pymol_process.cmd.reinitialize()


        # 清空 UI 相關欄位
        self.ui.textBrowser_workdirectory.clear()
        self.ui.tableWidget_show_receptor.setRowCount(0)
        self.ui.tableWidget_show_refligands.setRowCount(0)
        self.ui.tableWidget_show_ligands.setRowCount(0)
        self.ui.pushButton_uploadligands.setEnabled(False)
        
        
        
        self.ui.pushButton_setgridbox.setEnabled(False)
        self.ui.pushButton_setparameter.setEnabled(False)
        self.ui.pushButton_dockingbutton.setEnabled(False)
        
         