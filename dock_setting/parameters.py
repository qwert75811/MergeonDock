# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 16:14:37 2024

@author: Xhamrock Studio
"""

import os
from PyQt5.QtCore import QRegExp
from PyQt5.QtGui import QRegExpValidator




class Parameter_setting():
    def __init__(self, ui, pymol_process, all_parameters):
        self.ui = ui
        self.pymol_process = pymol_process  # 將 pymol_process 參數保存為類的屬性
        self.all_parameters = all_parameters

        self.lineEdit_validator()
    
    def set_gridbox_instance(self, gridbox_instance):
        self.gridbox_instance = gridbox_instance
        
    
    def lineEdit_validator(self):
        # 定义允许的字符（数字）
        regex = QRegExp("[0-9]*")
        # 创建 QRegExpValidator 验证器
        validator = QRegExpValidator(regex)
        # 为 QLineEdit 控件设置验证器
        self.ui.lineEdit_cpu.setValidator(validator)
        self.ui.lineEdit_seed.setValidator(validator)
    
    def switch_to_parameter_tab(self):
        # 解绑之前的连接
        self.ui.pushButton_setgridbox.clicked.disconnect()
        self.ui.pushButton_setparameter.clicked.disconnect()
    
        # 绑定 parameter 特定的功能
        self.ui.pushButton_setgridbox.clicked.connect(self.grid_button_save_button)
        self.ui.pushButton_setparameter.clicked.connect(self.parameter_button_cancel_button)
        # 更新按钮文本
        self.ui.pushButton_setgridbox.setText("Save")
        self.ui.pushButton_setparameter.setText("Cancel")
        # 切换到 parameter 页面
        self.ui.stackedWidget.setCurrentWidget(self.ui.page_parameters_setting)
        
        
        
    def grid_button_save_button(self):
        current_widget = self.ui.stackedWidget.currentWidget() 
        if current_widget == self.ui.page_parameters_setting:
            select_scoring_function = self.ui.buttonGroup_scoring_function.checkedButton()
            if select_scoring_function == self.ui.radioButton_vina:
                self.all_parameters.scoring_function = "vina"
            elif select_scoring_function == self.ui.radioButton_autodock4:
                self.all_parameters.scoring_function = "autodock4"
            
            self.all_parameters.exhaustiveness = self.ui.spinBox_exhaustiveness.value()
            self.all_parameters.poses = self.ui.spinBox_pose.value()
            
            select_verbosity = self.ui.buttonGroup_verbosity.checkedButton()
            if select_verbosity == self.ui.radioButton_verbosity_1:
                self.all_parameters.verbosity_value = 1
            elif select_verbosity == self.ui.radioButton_verbosity_2:
                self.all_parameters.verbosity_value = 2
            
            if self.ui.lineEdit_seed.text() == None or self.ui.lineEdit_seed.text() == "":
                self.all_parameters.seed_value = 0
            else:
                self.all_parameters.seed_value = self.ui.lineEdit_seed.text()
                
            if self.ui.lineEdit_cpu.text() == None or self.ui.lineEdit_cpu.text() == "":
                self.all_parameters.cpu_value = 0
            else:    
                self.all_parameters.cpu_value = self.ui.lineEdit_cpu.text()
                
            self.ui.stackedWidget.setCurrentWidget(self.ui.page_show_all_input)
            self.ui.pushButton_setgridbox.clicked.disconnect()
            self.ui.pushButton_setparameter.clicked.disconnect()
            

            self.ui.pushButton_setgridbox.clicked.connect(self.gridbox_instance.switch_to_gridbox_tab)
            self.ui.pushButton_setparameter.clicked.connect(self.switch_to_parameter_tab)

            self.ui.pushButton_setgridbox.setText("Grid box")
            self.ui.pushButton_setparameter.setText("Parameters")
            
            
        
    
    def parameter_button_cancel_button(self):
        current_widget = self.ui.stackedWidget.currentWidget()
        self.pymol_process.cmd.delete("gridbox")

        if current_widget == self.ui.page_parameters_setting:
            if self.all_parameters.scoring_function == "vina":
                self.ui.radioButton_vina.setChecked(True)
                self.ui.radioButton_autodock4.setChecked(False)
            elif self.all_parameters.scoring_function == "autodock4":
                self.ui.radioButton_vina.setChecked(False)
                self.ui.radioButton_autodock4.setChecked(True)
                
            self.ui.spinBox_exhaustiveness.setValue(self.all_parameters.exhaustiveness)
            self.ui.spinBox_pose.setValue(self.all_parameters.poses)
            
            if self.all_parameters.verbosity_value == 1:
                self.ui.radioButton_verbosity_1.setChecked(True)
                self.ui.radioButton_verbosity_2.setChecked(False)
            elif self.all_parameters.verbosity_value == 2:
                self.ui.radioButton_verbosity_1.setChecked(False)
                self.ui.radioButton_verbosity_2.setChecked(True)
            
            if self.all_parameters.cpu_value == 0:
                self.ui.lineEdit_cpu.setPlaceholderText("Auto")
            else:
                self.ui.lineEdit_cpu.setText(str(self.all_parameters.cpu_value))
                
            if self.all_parameters.seed_value == 0:
                self.ui.lineEdit_seed.setPlaceholderText("Random")
            else:
                self.ui.lineEdit_seed.setText(str(self.all_parameters.seed_value))
            
            
            
            self.ui.stackedWidget.setCurrentWidget(self.ui.page_show_all_input)
            self.ui.pushButton_setgridbox.clicked.disconnect()
            self.ui.pushButton_setparameter.clicked.disconnect()
            
            self.ui.pushButton_setgridbox.clicked.connect(self.gridbox_instance.switch_to_gridbox_tab)
            self.ui.pushButton_setparameter.clicked.connect(self.switch_to_parameter_tab)

            self.ui.pushButton_setgridbox.setText("Grid box")
            self.ui.pushButton_setparameter.setText("Parameters")
              
            
            
            
                
       
            
            
            
        