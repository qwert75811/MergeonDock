# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:19:20 2024

@author: Xhamrock Studio
"""

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QHeaderView, QTableWidgetItem, QMessageBox
from PyQt5.QtCore import Qt, QRegExp, QEvent
from PyQt5.QtGui import QRegExpValidator

import os
import subprocess
import re

from MergeonDock.menu import advance_setting_ui


class Menu_option_advance_setting():
    def __init__(self, ui, all_parameters):
        self.ui = ui
        self.ui.actionAdvance_setting.triggered.connect(self.open_own_window)
        self.all_parameters = all_parameters
        
        
    def open_own_window(self):
        self.own_window = Advance_setting_window(self.all_parameters)
        self.own_window.exec_()


class Advance_setting_window(QDialog):
    def __init__(self, all_parameters):
        super().__init__()
        self.ui_dialog = advance_setting_ui.Ui_Dialog_advance_setting()
        self.ui_dialog.setupUi(self)
        self.all_parameters = all_parameters
        self.previous_settings()
        self.ui_dialog.lineEdit_lig_ad4_I.installEventFilter(self)
        self.ad4_lig_I_groups = []
        self.ad4_lig_I_current_group = []
        self.ui_setting()
        
        
    def ui_setting(self):
        #Receptor preparation tab
        self.ui_dialog.radioButton_rec_ad4.toggled.connect(lambda: self.ui_dialog.stackedWidget_rec_ad4_opt_parameters.setCurrentIndex(0))
        self.ui_dialog.radioButton_rec_meeko.toggled.connect(lambda: self.ui_dialog.stackedWidget_rec_ad4_opt_parameters.setCurrentIndex(1))
        
        #Ligands preparation tab
        self.ui_dialog.radioButton_lig_ad4.toggled.connect(lambda: self.ui_dialog.stackedWidget_lig_ad4_opt_parameters.setCurrentIndex(0))
        self.ui_dialog.radioButton_lig_meeko.toggled.connect(lambda: self.ui_dialog.stackedWidget_lig_ad4_opt_parameters.setCurrentIndex(1))
        
        
        #Receptor optional window ad4
        self.ui_dialog.checkBox_rec_ad4_opt_parameters.toggled.connect(
            lambda checked: (
            self.ui_dialog.checkBox_rec_ad4_A.setEnabled(checked),
            self.ui_dialog.checkBox_rec_ad4_C.setEnabled(checked),
            self.ui_dialog.groupBox_rec_ad4_U_setting.setEnabled(checked),
            self.ui_dialog.checkBox_rec_ad4_d.setEnabled(checked),
            self.ui_dialog.checkBox_rec_ad4_e.setEnabled(checked),
            self.ui_dialog.checkBox_rec_ad4_p.setEnabled(checked),
            self.ui_dialog.checkBox_rec_ad4_w.setEnabled(checked)
            )
        )
        
        self.ui_dialog.checkBox_rec_ad4_A.toggled.connect(lambda checked:self.ui_dialog.comboBox_rec_ad4_A_setting.setEnabled(checked))
        self.ui_dialog.checkBox_rec_ad4_d.toggled.connect(lambda checked:self.ui_dialog.lineEdit_rec_ad4_d.setEnabled(checked))
        self.ui_dialog.checkBox_rec_ad4_p.toggled.connect(lambda checked:self.ui_dialog.lineEdit_rec_ad4_p.setEnabled(checked))
        
        
        
        #Receptor optional window meeko
        self.ui_dialog.checkBox_rec_meeko_opt_parameters.toggled.connect(
            lambda checked: (
            )
        )       
        
        #Ligands optional window ad4
        self.ui_dialog.checkBox_lig_ad4_opt_parameters.toggled.connect(
            lambda checked: (
            self.ui_dialog.checkBox_lig_ad4_A.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_C.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_Z.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_g.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_s.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_w.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_F.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_R.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_I.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_p.setEnabled(checked),
            self.ui_dialog.checkBox_lig_ad4_d.setEnabled(checked),
            self.ui_dialog.groupBox_lig_ad4_U_setting.setEnabled(checked),
            self.ui_dialog.groupBox_lig_ad4_B_setting.setEnabled(checked)
            )
        )
        
        self.ui_dialog.checkBox_lig_ad4_A.toggled.connect(lambda checked:self.ui_dialog.comboBox_lig_ad4_A_setting.setEnabled(checked))
        self.ui_dialog.checkBox_lig_ad4_R.toggled.connect(lambda checked:self.ui_dialog.lineEdit_lig_ad4_R.setEnabled(checked))
        self.ui_dialog.checkBox_lig_ad4_I.toggled.connect(lambda checked:self.ui_dialog.lineEdit_lig_ad4_I.setEnabled(checked))
        self.ui_dialog.checkBox_lig_ad4_p.toggled.connect(lambda checked:self.ui_dialog.lineEdit_lig_ad4_p.setEnabled(checked))
        self.ui_dialog.checkBox_lig_ad4_d.toggled.connect(lambda checked:self.ui_dialog.lineEdit_lig_ad4_d.setEnabled(checked))
        
        
        #Ligands optional window meeko
        self.ui_dialog.checkBox_meeko_lig_opt_parameters.toggled.connect(
            lambda checked: (
            )
        )     
   
        self.ad4_p_lineEdit_validator()
        self.ad4_lig_R_lineEdit_validator()
        self.ui_dialog.pushButton_Cancel.clicked.connect(self.cancel_button)
        self.ui_dialog.pushButton_Save.clicked.connect(self.save_button)
        
        
    
    def eventFilter(self, source, event):                                                              #這段的邏輯有一個問題是不支援左右鍵移動做修改
        if event.type() == QEvent.KeyPress and source is self.ui_dialog.lineEdit_lig_ad4_I:
            valid_characters = "0123456789,"
            
            char = event.text()
            
            if char and char not in valid_characters and event.key() not in (Qt.Key_Backspace, Qt.Key_Delete, Qt.Key_Left, Qt.Key_Right):
                return True  # Ignore invalid characters
        
            if char.isdigit():
                digit = event.text()  # 获取按下的键的字符 
                self.ad4_lig_I_current_group.append(digit)  #將當前輸入數字放入目前列表
                self.update_display()
                return True
            elif event.key() == Qt.Key_Comma:
                if self.ad4_lig_I_current_group:
                    self.ad4_lig_I_groups.append(''.join(self.ad4_lig_I_current_group)) #將目前列表整合成一組,並放到另一個列表中代表一組數字
                    self.ad4_lig_I_current_group = [] #清空目前列表,等待下一組數字輸入
                self.update_display()
                return True
            elif event.key() == Qt.Key_Backspace:
                # 捕获退格键事件
                if self.ad4_lig_I_current_group:
                    self.ad4_lig_I_current_group.pop()
                elif self.ad4_lig_I_groups:
                    self.ad4_lig_I_current_group = list(self.ad4_lig_I_groups.pop())
                self.update_display()
                return True
        return super(Advance_setting_window, self).eventFilter(source, event)

    def update_display(self):
        self.ad4_lig_I_formatted_groups = []
        
        # Format groups already entered
        for i in range(0, len(self.ad4_lig_I_groups), 2):        #已經轉化到self.group的部分,倆倆一組放到formatted group做顯示
            if i+1 < len(self.ad4_lig_I_groups):
                self.ad4_lig_I_formatted_groups.append(f'({self.ad4_lig_I_groups[i]},{self.ad4_lig_I_groups[i+1]})')   #數量超過2組就倆倆顯示
            else:
                self.ad4_lig_I_formatted_groups.append(f'({self.ad4_lig_I_groups[i]})')                      #數量不足2就只先顯示一組
    
        # Format current group being entered          
        if self.ad4_lig_I_current_group:                                                      #尚未用逗點觸發使current group轉化到self.group, 如果current group存在則繼續顯示
            if len(self.ad4_lig_I_groups) % 2 == 0:                                           #先檢查self.group是否倆倆一組,若是且current group存在則代表有新的一組正在輸入的數字,加入到formatted group一起顯示
                self.ad4_lig_I_formatted_groups.append(f'({"".join(self.ad4_lig_I_current_group)})')
            else:
                if len(self.ad4_lig_I_current_group) > 1:                                                                                           #正在輸入尚未轉化進self.group的部分,檢查是否長度>1,若有代表有正在輸入的數字
                    self.ad4_lig_I_formatted_groups[-1] = f'{self.ad4_lig_I_formatted_groups[-1][:-1]},{"".join(self.ad4_lig_I_current_group)})'    #把尚未進入formatted group的數字放入,把formatted group最後一項的新增上current group正在輸入的內容
                else:                                                                                                                               #[:-1] 刪掉最後一個括號去接上current group
                    self.ad4_lig_I_formatted_groups[-1] = f'{self.ad4_lig_I_formatted_groups[-1][:-1]},{self.ad4_lig_I_current_group[0]})'
 
        
        self.ui_dialog.lineEdit_lig_ad4_I.setText(','.join(self.ad4_lig_I_formatted_groups))
        
        




    
    
    def previous_settings(self):
        #記憶保存之前的設定
        #-----------------autodock4_receptor----------------------------------------------------------------------
        if self.all_parameters.receptor_prepare_method == "ad4":
            self.ui_dialog.radioButton_rec_ad4.setChecked(True)
        if self.all_parameters.receptor_prepare_opt_switch == True:
            self.ui_dialog.checkBox_rec_ad4_opt_parameters.setChecked(True)
            self.ui_dialog.checkBox_rec_ad4_A.setEnabled(True)
            self.ui_dialog.checkBox_rec_ad4_C.setEnabled(True)
            self.ui_dialog.groupBox_rec_ad4_U_setting.setEnabled(True)
            self.ui_dialog.checkBox_rec_ad4_d.setEnabled(True)
            self.ui_dialog.checkBox_rec_ad4_e.setEnabled(True)
            self.ui_dialog.checkBox_rec_ad4_p.setEnabled(True)
            self.ui_dialog.checkBox_rec_ad4_w.setEnabled(True)
            
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] != "":
            self.ui_dialog.checkBox_rec_ad4_A.setChecked(True)
            self.ui_dialog.comboBox_rec_ad4_A_setting.setEnabled(True)
            self.ui_dialog.comboBox_rec_ad4_A_setting.setCurrentIndex(self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A_combobox"])
        
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_C"] != "":
            self.ui_dialog.checkBox_rec_ad4_C.setChecked(True)
        
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_e"] != "":
            self.ui_dialog.checkBox_rec_ad4_e.setChecked(True)
        
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_w"] != "":
            self.ui_dialog.checkBox_rec_ad4_w.setChecked(True)
        
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p"] != "":
            self.ui_dialog.checkBox_rec_ad4_p.setChecked(True)
            self.ui_dialog.lineEdit_rec_ad4_p.setEnabled(True)
            self.ui_dialog.lineEdit_rec_ad4_p.setText(self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p_lineedit"])
        
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d"] != "":
            self.ui_dialog.checkBox_rec_ad4_d.setChecked(True)
            self.ui_dialog.lineEdit_rec_ad4_d.setEnabled(True)
            self.ui_dialog.lineEdit_rec_ad4_d.setText(self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d_lineedit"])
        
        if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U"] != "":
            self.ui_dialog.groupBox_rec_ad4_U_setting.setChecked(True)
            if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_nphs"] != "":
                self.ui_dialog.checkBox_rec_ad4_U_nphs.setChecked(True)
            else:
                self.ui_dialog.checkBox_rec_ad4_U_nphs.setChecked(False)
            if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_lps"] != "":
                self.ui_dialog.checkBox_rec_ad4_U_lps.setChecked(True)
            else:
                self.ui_dialog.checkBox_rec_ad4_U_lps.setChecked(False)
            if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_waters"] != "":
                self.ui_dialog.checkBox_rec_ad4_U_waters.setChecked(True)
            else:
                self.ui_dialog.checkBox_rec_ad4_U_waters.setChecked(False)
            if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_nonstdres"] != "":
                self.ui_dialog.checkBox_rec_ad4_U_nonstdres.setChecked(True)
            else:
                self.ui_dialog.checkBox_rec_ad4_U_nonstdres.setChecked(False)
        
         #----------------autodock4_ligands----------------------------------------------------------------------
        if self.all_parameters.ligands_prepare_method == "ad4":
            self.ui_dialog.radioButton_lig_ad4.setChecked(True)
        if self.all_parameters.ligands_prepare_opt_switch == True:
            self.ui_dialog.checkBox_lig_ad4_opt_parameters.setChecked(True)
            self.ui_dialog.checkBox_lig_ad4_A.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_C.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_Z.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_g.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_s.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_w.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_F.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_R.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_I.setEnabled(True) 
            self.ui_dialog.checkBox_lig_ad4_d.setEnabled(True)
            self.ui_dialog.checkBox_lig_ad4_p.setEnabled(True)
            self.ui_dialog.groupBox_lig_ad4_U_setting.setEnabled(True)
            self.ui_dialog.groupBox_lig_ad4_B_setting.setEnabled(True)
            
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A"] != "":
            self.ui_dialog.checkBox_lig_ad4_A.setChecked(True)
            self.ui_dialog.comboBox_lig_ad4_A_setting.setEnabled(True)
            self.ui_dialog.comboBox_lig_ad4_A_setting.setCurrentIndex(self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A_combobox"])
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_C"] != "":
            self.ui_dialog.checkBox_lig_ad4_C.setChecked(True)
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_Z"] != "":
            self.ui_dialog.checkBox_lig_ad4_Z.setChecked(True)
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_g"] != "":
            self.ui_dialog.checkBox_lig_ad4_g.setChecked(True)
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_s"] != "":
            self.ui_dialog.checkBox_lig_ad4_s.setChecked(True)
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_w"] != "":
            self.ui_dialog.checkBox_lig_ad4_w.setChecked(True)
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_F"] != "":
            self.ui_dialog.checkBox_lig_ad4_F.setChecked(True)
        
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_R"] != "":
            self.ui_dialog.checkBox_lig_ad4_R.setChecked(True)
            self.ui_dialog.lineEdit_lig_ad4_R.setEnabled(True)
            self.ui_dialog.lineEdit_lig_ad4_R.setText(self.all_parameters.ligands_opt_parameters_dict["lig_ad4_R_lineedit"])
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_I"] != "":
            self.ui_dialog.checkBox_lig_ad4_I.setChecked(True)
            self.ui_dialog.lineEdit_lig_ad4_I.setEnabled(True)
            self.ui_dialog.lineEdit_lig_ad4_I.setText(self.all_parameters.ligands_opt_parameters_dict["lig_ad4_I_lineedit"])
            
        
        
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p"] != "":
            self.ui_dialog.checkBox_lig_ad4_p.setChecked(True)
            self.ui_dialog.lineEdit_lig_ad4_p.setEnabled(True)
            self.ui_dialog.lineEdit_lig_ad4_p.setText(self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p_lineedit"])
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d"] != "":
            self.ui_dialog.checkBox_lig_ad4_d.setChecked(True)
            self.ui_dialog.lineEdit_lig_ad4_d.setEnabled(True)
            self.ui_dialog.lineEdit_lig_ad4_d.setText(self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d_lineedit"])
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U"] != "":
            self.ui_dialog.groupBox_lig_ad4_U_setting.setChecked(True)
            if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U_nphs"] != "":
                self.ui_dialog.checkBox_lig_ad4_U_nphs.setChecked(True)
            else:
                self.ui_dialog.checkBox_lig_ad4_U_nphs.setChecked(False)
                
            if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U_lps"] != "":
                self.ui_dialog.checkBox_lig_ad4_U_lps.setChecked(True)
            else:
                self.ui_dialog.checkBox_lig_ad4_U_lps.setChecked(False)
                
        
        if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B"] != "":
            self.ui_dialog.groupBox_lig_ad4_B_setting.setChecked(True)
            if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_backbone"] != "":
                self.ui_dialog.checkBox_lig_ad4_B_backbone.setChecked(True)
            else:
                self.ui_dialog.checkBox_lig_ad4_B_backbone.setChecked(False)
                
            if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_amide"] != "":
                self.ui_dialog.checkBox_lig_ad4_B_amide.setChecked(True)
            else:
                self.ui_dialog.checkBox_lig_ad4_B_amide.setChecked(False)
            
            if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_guanidinium"] != "":
                self.ui_dialog.checkBox_lig_ad4_B_guanidinium.setChecked(True)
            else:
                self.ui_dialog.checkBox_lig_ad4_B_guanidinium.setChecked(False)
        
    
    
    def ad4_p_lineEdit_validator(self):
        # ad4 -p參數的輸入格式限制
        # 定义允许的字符（英文、数字、逗號、空白）
        regex = QRegExp("[a-zA-Z0-9, ]*")
        # 创建 QRegExpValidator 验证器
        validator = QRegExpValidator(regex, self)
        # 为 QLineEdit 控件设置验证器
        self.ui_dialog.lineEdit_rec_ad4_p.setValidator(validator)
        self.ui_dialog.lineEdit_lig_ad4_p.setValidator(validator)
    
    def ad4_lig_R_lineEdit_validator(self):
        # ad4 lig -R參數的輸入格式限制
        # 定义允许的字符（数字）
        regex = QRegExp("[0-9]*")
        # 创建 QRegExpValidator 验证器
        validator = QRegExpValidator(regex, self)
        # 为 QLineEdit 控件设置验证器
        self.ui_dialog.lineEdit_lig_ad4_R.setValidator(validator)
        
    
    
    
    
    def save_button(self):
        invaild_signal = False    
    #---------------------------------------receptor選ad4 方法-------------------------------------------------------------------------------
        if self.ui_dialog.radioButton_rec_ad4.isChecked():
            self.all_parameters.receptor_prepare_method = "ad4"
            #Optional 有沒有開
            if self.ui_dialog.checkBox_rec_ad4_opt_parameters.isChecked():
                self.all_parameters.receptor_prepare_opt_switch = True    
            else:
                self.all_parameters.receptor_prepare_opt_switch = False
                
            #receptor optional -A    
            if self.ui_dialog.checkBox_rec_ad4_A.isChecked():
                rec_ad4_A_opt_index = self.ui_dialog.comboBox_rec_ad4_A_setting.currentIndex()
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A_combobox"] = rec_ad4_A_opt_index
                if rec_ad4_A_opt_index == 0:
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] = "-A None"
                elif rec_ad4_A_opt_index == 1:
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] = "-A hydrogens"
                elif rec_ad4_A_opt_index == 2:
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] = "-A checkhydrogens"
                elif rec_ad4_A_opt_index == 3:
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] = "-A bonds_hydrogens"
                elif rec_ad4_A_opt_index == 4:
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] = "-A bonds"
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_A"] = ""
                
            #receptor optional -C  
            if self.ui_dialog.checkBox_rec_ad4_C.isChecked():
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_C"] = "-C"
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_C"] = ""
                
            #receptor optional -e
            if self.ui_dialog.checkBox_rec_ad4_e.isChecked():
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_e"] = "-e"
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_e"] = ""
                
            #receptor optional -w
            if self.ui_dialog.checkBox_rec_ad4_w.isChecked():
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_w"] = "-w"
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_w"] = ""
                
            #receptor optional -p
            if self.ui_dialog.checkBox_rec_ad4_p.isChecked():
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p_lineedit"] = self.ui_dialog.lineEdit_rec_ad4_p.text().strip()
                rec_ad4_p_input_format = re.compile(r"^[a-zA-Z0-9, ]+$") #手動再次檢查輸入格式
                if self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p_lineedit"] == "":
                    warning_window = QMessageBox()
                    warning_window.setIcon(QMessageBox.Warning)
                    warning_window.setWindowTitle("Warning")
                    warning_window.setText("'-p' Please fill in the content. eg: Zn,Fe")
                    warning_window.setStandardButtons(QMessageBox.Ok)
                    warning_window.exec_()
                    invaild_signal = True
                elif not rec_ad4_p_input_format.match(self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p_lineedit"]):
                    warning_window = QMessageBox()
                    warning_window.setIcon(QMessageBox.Warning)
                    warning_window.setWindowTitle("Warning: Invalid Input")
                    warning_window.setText("'-p' Input contains invalid characters. eg: Zn,Fe")
                    warning_window.setStandardButtons(QMessageBox.Ok)
                    warning_window.exec_()
                    invaild_signal = True
                else:
                    atom_command = []
                    rec_ad4_p_final_command = ""
                    for atom in self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p_lineedit"].split(","):
                        if atom != "":
                            atom_command.append(f"-p {atom}")
                    for command in atom_command:
                        rec_ad4_p_final_command += command + " "
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p"] = rec_ad4_p_final_command
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_p"] = ""
                      
            #receptor optional -d 
            if self.ui_dialog.checkBox_rec_ad4_d.isChecked():
                rec_ad4_d_input_text = self.ui_dialog.lineEdit_rec_ad4_d.text().strip()
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d_lineedit"] = rec_ad4_d_input_text
                
                if not rec_ad4_d_input_text:
                    warning_window = QMessageBox()
                    warning_window.setIcon(QMessageBox.Warning)
                    warning_window.setWindowTitle("Warning")
                    warning_window.setText(r"'-d' Please fill in the content. eg:\path\File name")
                    warning_window.setStandardButtons(QMessageBox.Ok)
                    warning_window.exec_()
                    invaild_signal = True
                else:
                    directory_path = os.path.normpath(os.path.dirname(rec_ad4_d_input_text))
                    
                    if directory_path == ".":
                        full_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, rec_ad4_d_input_text))
                        dot_index = rec_ad4_d_input_text.rfind(".")
                        if dot_index != 0 and dot_index != -1: 
                            self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d"] = f"-d {full_path}"
                        else:
                            self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d"] = f"-d {full_path}.txt"  
                    else:
                        if os.path.isdir(directory_path):
                            full_path = rec_ad4_d_input_text
                            dot_index = full_path.rfind(".")
                            if dot_index != 0 and dot_index != -1:
                                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d"] = f"-d {full_path}"
                            else:
                                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d"] = f"-d {full_path}.txt"

                        else:
                            warning_window = QMessageBox()
                            warning_window.setIcon(QMessageBox.Warning)
                            warning_window.setWindowTitle("Warning")
                            warning_window.setText(r"'-d' Directory is not exist. Please enter the correct path")
                            warning_window.setStandardButtons(QMessageBox.Ok)
                            warning_window.exec_()
                            invaild_signal = True
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_d"] = ""
                
                
            #receptor optional -U
            if self.ui_dialog.groupBox_rec_ad4_U_setting.isChecked():
                U_command = ["", "", "", ""]
                if self.ui_dialog.checkBox_rec_ad4_U_nphs.isChecked():
                    U_command[0]= "nphs"
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_nphs"] = "checked"
                else:
                    U_command[0]= ""
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_nphs"] = ""
                    
                if self.ui_dialog.checkBox_rec_ad4_U_lps.isChecked():
                    U_command[1]= "lps"
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_lphs"] = "checked"
                else:
                    U_command[1]= ""
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_lphs"] = ""
                
                if self.ui_dialog.checkBox_rec_ad4_U_waters.isChecked():
                    U_command[2]= "waters"
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_waters"] = "checked"
                else:
                    U_command[2]= ""
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_waters"] = ""
                
                if self.ui_dialog.checkBox_rec_ad4_U_nonstdres.isChecked():
                    U_command[3]= "nonstdres"
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_nonstdres"] = "checked"
                else:
                    U_command[3]= ""
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U_nonstdres"] = ""
                
                U_final_command = ""
                if U_command[0] != "":
                    U_final_command += U_command[0] + "_"
                if U_command[1] != "":
                    U_final_command += U_command[1] + "_"
                if U_command[2] != "":
                    U_final_command += U_command[2] + "_"
                if U_command[3] != "":
                    U_final_command += U_command[3] + "_"
                    
                if U_final_command == "":
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U"] = ""
                else:
                    U_final_command = U_final_command.rstrip("_")  # 移除尾部多余的下划线
                    self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U"] = f"-U {U_final_command}"
            else:
                self.all_parameters.receptor_opt_parameters_dict["rec_ad4_U"] = ""
                
            #將存取的自訂參數合併成指令
            prepare_receptor_custom_command = []
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_A'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_A'])
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_C'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_C'])
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_e'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_e'])
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_w'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_w'])
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_p'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_p'])
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_d'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_d'])
            if self.all_parameters.receptor_opt_parameters_dict['rec_ad4_U'] != "":
                prepare_receptor_custom_command.append(self.all_parameters.receptor_opt_parameters_dict['rec_ad4_U'])
            
            self.all_parameters.autodock_prepare_receptor_custom_command = " ".join(prepare_receptor_custom_command)
            
        
            
    #-----------------------------------------receptor選 meeko 方法------------------------------------------------------------------------------------------
        if self.ui_dialog.radioButton_rec_meeko.isChecked():
            self.all_parameters.receptor_prepare_method = "meeko"
        
        else:
            self.all_parameters.receptor_prepare_method = "ad4"
            
    
    #-----------------------------------------ligands選 ad4 方法---------------------------------------------------------------------------------------------
        if self.ui_dialog.radioButton_lig_ad4.isChecked():
            self.all_parameters.ligands_prepare_method = "ad4"
            #Optional 有沒有開
            if self.ui_dialog.checkBox_lig_ad4_opt_parameters.isChecked():
                self.all_parameters.ligands_prepare_opt_switch = True    
            else:
                self.all_parameters.ligands_prepare_opt_switch = False
                
        #ligands optional -A    
        if self.ui_dialog.checkBox_lig_ad4_A.isChecked():
            lig_ad4_A_opt_index = self.ui_dialog.comboBox_lig_ad4_A_setting.currentIndex()
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A_combobox"] = lig_ad4_A_opt_index
            if lig_ad4_A_opt_index == 0:
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A"] = "-A None"
            elif lig_ad4_A_opt_index == 1:
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A"] = "-A hydrogens"
            elif lig_ad4_A_opt_index == 2:
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A"] = "-A bonds_hydrogens"
            elif lig_ad4_A_opt_index == 3:
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A"] = "-A bonds"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_A"] = ""
            
        #ligands optional -C  
        if self.ui_dialog.checkBox_lig_ad4_C.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_C"] = "-C"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_C"] = ""
            
        #ligands optional -Z
        if self.ui_dialog.checkBox_lig_ad4_Z.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_Z"] = "-Z"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_Z"] = ""
        
        #ligands optional -g
        if self.ui_dialog.checkBox_lig_ad4_g.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_g"] = "-g"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_g"] = ""
        
        #ligands optional -s
        if self.ui_dialog.checkBox_lig_ad4_s.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_s"] = "-s"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_s"] = ""
            
        #ligands optional -w
        if self.ui_dialog.checkBox_lig_ad4_w.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_w"] = "-w"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_w"] = ""
        
        #ligands optional -F
        if self.ui_dialog.checkBox_lig_ad4_F.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_F"] = "-F"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_F"] = ""
            
        #ligands optional -R
        if self.ui_dialog.checkBox_lig_ad4_R.isChecked():
            lig_ad4_R_input_text = self.ui_dialog.lineEdit_lig_ad4_R.text().strip()
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_R_lineedit"] = lig_ad4_R_input_text
            
            if not lig_ad4_R_input_text:
                warning_window = QMessageBox()
                warning_window.setIcon(QMessageBox.Warning)
                warning_window.setWindowTitle("Warning")
                warning_window.setText(r"'-R' Please fill in the content. eg:15 (Number only)")
                warning_window.setStandardButtons(QMessageBox.Ok)
                warning_window.exec_()
                invaild_signal = True
            else:
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_R"] = f"-R {lig_ad4_R_input_text}"

        #ligands optional -I
        if self.ui_dialog.checkBox_lig_ad4_I.isChecked():
            lig_ad4_I_input_text = self.ui_dialog.lineEdit_lig_ad4_I.text().strip()
            self.lig_ad4_I_lineedit_savedtext = ','.join(self.ad4_lig_I_formatted_groups)
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_I_lineedit"] = self.lig_ad4_I_lineedit_savedtext
            
            if not lig_ad4_I_input_text:
                warning_window = QMessageBox()
                warning_window.setIcon(QMessageBox.Warning)
                warning_window.setWindowTitle("Warning")
                warning_window.setText(r"'-I' Please fill in the content. eg:15,25 (Number, comma only)")
                warning_window.setStandardButtons(QMessageBox.Ok)
                warning_window.exec_()
                invaild_signal = True
            else:
                total_atom_list = []
                for atom_text in self.ad4_lig_I_formatted_groups:
                    number = atom_text.strip("()").split(",") 
                    for num in number:
                        total_atom_list.append(int(num))

                if len(total_atom_list) %2 != 0:
                    warning_window = QMessageBox()
                    warning_window.setIcon(QMessageBox.Warning)
                    warning_window.setWindowTitle("Warning")
                    warning_window.setText(r"'-I' Atom number must be in pairs. eg:(15,25),(120,234)")
                    warning_window.setStandardButtons(QMessageBox.Ok)
                    warning_window.exec_()
                    invaild_signal = True    
                else:
                    atom_pairs_list = []
                    for atom_pair in self.ad4_lig_I_formatted_groups:
                        atom_number = atom_pair.strip("()").split(",")
                        atom_pair_num = [int(num) for num in atom_number]
                        atom_pairs_list.append(tuple(atom_pair_num))
                    
                    
                    atom_command_format = ""
                    for tuple_pair in atom_pairs_list:
                        atom_command_format += f"{tuple_pair[0]}_{tuple_pair[1]}_"
                        
                    atom_command_format = atom_command_format.rstrip("_")
                    self.all_parameters.ligands_opt_parameters_dict["lig_ad4_I"] = f"-I {atom_command_format}"
                    
                    
                    

        
        #ligands optional -p
        if self.ui_dialog.checkBox_lig_ad4_p.isChecked():
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p_lineedit"] = self.ui_dialog.lineEdit_lig_ad4_p.text().strip()
            lig_ad4_p_input_format = re.compile(r"^[a-zA-Z0-9, ]+$") #手動再次檢查輸入格式
            if self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p_lineedit"] == "":
                warning_window = QMessageBox()
                warning_window.setIcon(QMessageBox.Warning)
                warning_window.setWindowTitle("Warning")
                warning_window.setText("'-p' Please fill in the content. eg: Zn,Fe")
                warning_window.setStandardButtons(QMessageBox.Ok)
                warning_window.exec_()
                invaild_signal = True
            elif not lig_ad4_p_input_format.match(self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p_lineedit"]):
                warning_window = QMessageBox()
                warning_window.setIcon(QMessageBox.Warning)
                warning_window.setWindowTitle("Warning: Invalid Input")
                warning_window.setText("'-p' Input contains invalid characters. eg: Zn,Fe")
                warning_window.setStandardButtons(QMessageBox.Ok)
                warning_window.exec_()
                invaild_signal = True
            else:
                atom_command = []
                lig_ad4_p_final_command = ""
                for atom in self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p_lineedit"].split(","):
                    if atom != "":
                        atom_command.append(f"-p {atom}")
                for command in atom_command:
                    lig_ad4_p_final_command += command + " "
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p"] = lig_ad4_p_final_command
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_p"] = ""
                  
        #ligands optional -d 
        if self.ui_dialog.checkBox_lig_ad4_d.isChecked():
            lig_ad4_d_input_text = self.ui_dialog.lineEdit_lig_ad4_d.text().strip()
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d_lineedit"] = lig_ad4_d_input_text
            
            if not lig_ad4_d_input_text:
                warning_window = QMessageBox()
                warning_window.setIcon(QMessageBox.Warning)
                warning_window.setWindowTitle("Warning")
                warning_window.setText(r"'-d' Please fill in the content. eg:\path\File name")
                warning_window.setStandardButtons(QMessageBox.Ok)
                warning_window.exec_()
                invaild_signal = True
            else:
                directory_path = os.path.normpath(os.path.dirname(lig_ad4_d_input_text))
                
                if directory_path == ".":
                    full_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, lig_ad4_d_input_text))
                    dot_index = lig_ad4_d_input_text.rfind(".")
                    if dot_index != 0 and dot_index != -1: 
                        self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d"] = f"-d {full_path}"
                    else:
                        self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d"] = f"-d {full_path}.txt"  
                else:
                    if os.path.isdir(directory_path):
                        full_path = lig_ad4_d_input_text
                        dot_index = full_path.rfind(".")
                        if dot_index != 0 and dot_index != -1:
                            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d"] = f"-d {full_path}"
                        else:
                            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d"] = f"-d {full_path}.txt"

                    else:
                        warning_window = QMessageBox()
                        warning_window.setIcon(QMessageBox.Warning)
                        warning_window.setWindowTitle("Warning")
                        warning_window.setText(r"'-d' Directory is not exist. Please enter the correct path")
                        warning_window.setStandardButtons(QMessageBox.Ok)
                        warning_window.exec_()
                        invaild_signal = True
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_d"] = ""
        
        #ligands optional -U
        if self.ui_dialog.groupBox_lig_ad4_U_setting.isChecked():
            U_command = ["", ""]
            if self.ui_dialog.checkBox_lig_ad4_U_nphs.isChecked():
                U_command[0]= "nphs"
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U_nphs"] = "checked"
            else:
                U_command[0]= ""
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U_nphs"] = ""
                
            if self.ui_dialog.checkBox_lig_ad4_U_lps.isChecked():
                U_command[1]= "lps"
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U_lps"] = "checked"
            else:
                U_command[1]= ""
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U_lps"] = ""
            
            
            
            U_final_command = ""
            if U_command[0] != "":
                U_final_command += U_command[0] + "_"
            if U_command[1] != "":
                U_final_command += U_command[1] + "_"
           
                
            if U_final_command == "":
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U"] = ""
            else:
                U_final_command = U_final_command.rstrip("_")  # 移除尾部多余的下划线
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U"] = f"-U {U_final_command}"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_U"] = ""
        
        
        #ligands optional -B
        if self.ui_dialog.groupBox_lig_ad4_B_setting.isChecked():
            lig_ad4_B_command = ["", "", ""]
            if self.ui_dialog.checkBox_lig_ad4_B_backbone.isChecked():
                lig_ad4_B_command[0]= "backbone"
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_backbone"] = "checked"
            else:
                lig_ad4_B_command[0]= ""
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_backbone"] = ""
                
            if self.ui_dialog.checkBox_lig_ad4_B_amide.isChecked():
                lig_ad4_B_command[1]= "amide"
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_amide"] = "checked"
            else:
                lig_ad4_B_command[1]= ""
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_amide"] = ""
            
            if self.ui_dialog.checkBox_lig_ad4_B_guanidinium.isChecked():
                lig_ad4_B_command[2]= "guanidinium"
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_guanidinium"] = "checked"
            else:
                lig_ad4_B_command[2]= ""
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B_guanidinium"] = ""
            
            
            
            lig_ad4_B_final_command = ""
            if lig_ad4_B_command[0] != "":
                lig_ad4_B_final_command += lig_ad4_B_command[0] + " "
            if lig_ad4_B_command[1] != "":
                lig_ad4_B_final_command += lig_ad4_B_command[1] + " "
            if lig_ad4_B_command[2] != "":
                lig_ad4_B_final_command += lig_ad4_B_command[2] + " "
           
                
            if lig_ad4_B_final_command == "":
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B"] = ""
            else:
                lig_ad4_B_final_command = lig_ad4_B_final_command.rstrip(" ")  # 移除尾部多余的下划线
                self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B"] = f"-B {lig_ad4_B_final_command}"
        else:
            self.all_parameters.ligands_opt_parameters_dict["lig_ad4_B"] = ""
        
        
        
        #將存取的自訂參數合併成指令
        prepare_ligands_custom_command = []
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_A'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_A'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_C'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_C'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_Z'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_Z'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_g'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_g'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_s'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_s'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_w'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_w'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_F'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_F'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_p'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_p'])  
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_d'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_d'])
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_R'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_R'])
            
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_I'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_I'])
        
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_U'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_U'])
        
        if self.all_parameters.ligands_opt_parameters_dict['lig_ad4_B'] != "":
            prepare_ligands_custom_command.append(self.all_parameters.ligands_opt_parameters_dict['lig_ad4_B'])
        
        self.all_parameters.autodock_prepare_ligands_custom_command = " ".join(prepare_ligands_custom_command)
        
        print(self.all_parameters.autodock_prepare_ligands_custom_command)
        
        if invaild_signal == False:
            self.close()
            
            
    def cancel_button(self):
        self.close()
    