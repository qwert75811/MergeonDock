# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:10:36 2024

@author: Xhamrock Studio
"""

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog, QHeaderView, QTableWidgetItem, QMenu
from PyQt5.QtCore import Qt

from openbabel import openbabel

import os
import sys
from MergeonDock.menu import File_format_converter_ui
import subprocess

class Menu_option_file_convert():
    def __init__(self, ui, all_parameters):
        self.ui = ui
        self.all_parameters = all_parameters
        self.ui.actionFile_format_convert_OpenBabel_3_1.triggered.connect(self.open_own_window)
        
    def open_own_window(self):
        self.own_window = Format_convert_window(self.all_parameters)
        self.own_window.exec_()


class Format_convert_window(QDialog):
    def __init__(self, all_parameters):
        super().__init__()
        self.ui_dialog = File_format_converter_ui.Ui_Dialog()
        self.ui_dialog.setupUi(self)
        self.all_parameters = all_parameters
        
        self.files_input_detail = {}
        
        self.ui_setting()
        
        
        #觸發按鈕
        self.ui_dialog.pushButton_file_upload.clicked.connect(self.files_input)
        self.ui_dialog.pushButton_convert.clicked.connect(self.files_convert)
        self.ui_dialog.pushButton_clearall.clicked.connect(self.clear_all)
        
        #設置右鍵菜單
        self.ui_dialog.tableWidget_file_list.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui_dialog.tableWidget_file_list.customContextMenuRequested.connect(self.right_click_menu)
        
        self.ui_dialog.tableWidget_file_list.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)      #指定用戶在點擊單元格時應選擇整行
        self.ui_dialog.tableWidget_file_list.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        
        
    def ui_setting(self):
        #設定成水平伸展貼合table widget
        self.ui_dialog.tableWidget_file_list.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        header_show_extention =  self.ui_dialog.tableWidget_file_list.horizontalHeader()        
        header_show_extention.setSectionResizeMode(1, QHeaderView.ResizeToContents)      # 第1列根據內容調整
        
        # 设置表格的单元格不可编辑
        self.ui_dialog.tableWidget_file_list.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        # 調整 windowFlags 屬性
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowMinMaxButtonsHint | QtCore.Qt.WindowCloseButtonHint)
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowContextHelpButtonHint)
        
        if len(self.files_input_detail) > 0:
            self.ui_dialog.pushButton_convert.setEnabled(True)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(255, 255, 255); background-color: rgb(60, 60, 60)")
        else:
            self.ui_dialog.pushButton_convert.setEnabled(False)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(100, 100, 100); background-color: rgb(30, 30, 30)")
        
        
    def files_input(self):
        files_input_path_raw = QtWidgets.QFileDialog.getOpenFileNames()      
        files_path = files_input_path_raw[0]
   
        file_amount = len(self.files_input_detail)
 
        for path in files_path:
            file_basename = os.path.basename(path)
            if file_basename in self.files_input_detail:
                continue
            
            file_basename = os.path.basename(path)
            file_name, file_extention = os.path.splitext(file_basename)
            file_extention = file_extention[1:]
            self.files_input_detail[file_basename] = {"name":file_name,"path":path, "extention":file_extention}
            
            self.ui_dialog.tableWidget_file_list.setRowCount(file_amount + 1)
            self.ui_dialog.tableWidget_file_list.setItem(file_amount, 0, QTableWidgetItem(file_name))  
            self.ui_dialog.tableWidget_file_list.setItem(file_amount, 1, QTableWidgetItem(file_extention))  
            file_amount += 1
        
        if len(self.files_input_detail) > 0:
            self.ui_dialog.pushButton_convert.setEnabled(True)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(255, 255, 255); background-color: rgb(60, 60, 60)")
        else:
            self.ui_dialog.pushButton_convert.setEnabled(False)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(100, 100, 100); background-color: rgb(30, 30, 30)")
            
 
        
    def files_convert(self):
        output_format_select = self.ui_dialog.comboBox_format.currentText()
    
        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetOutFormat(str(output_format_select))
    
        error_file_message = []
    
        # 检查用户是否选择了自定义输出文件夹
        if self.ui_dialog.checkBox_customize_output_folder.checkState() == Qt.Checked:
            # 获取用户选择的目录路径
            directory_path = QtWidgets.QFileDialog.getExistingDirectory()
        else:
            # 使用默认的输入文件所在目录
            directory_path = None
    
        # 遍历文件并进行转换
        for file_basename in self.files_input_detail:
            file_name = self.files_input_detail[file_basename]["name"]
            file_path = self.files_input_detail[file_basename]["path"]
            file_extention = self.files_input_detail[file_basename]["extention"]
    
            if directory_path:
                output_dir = directory_path
            else:
                output_dir = os.path.normpath(os.path.dirname(file_path))
    
            # 输出文件的路径
            output_path = os.path.normpath(os.path.join(output_dir, f"{file_name}.{str(output_format_select)}"))
    
            # 尝试进行文件转换
            try:
                ob_conversion.SetInFormat(str(file_extention))
                molecular = openbabel.OBMol()
                if not ob_conversion.ReadFile(molecular, file_path):
                    raise Exception("Failed to read input file")
                if not ob_conversion.WriteFile(molecular, output_path):
                    raise Exception("Failed to write output file")
            except Exception as e:
                error_file_message.append(file_path)
    
        # 如果有转换失败的文件，显示错误对话框
        if error_file_message:
            error_message = "Failed to convert the following files:\n" + "\n".join(error_file_message)
            error_window = QtWidgets.QMessageBox()
            error_window.setIcon(QtWidgets.QMessageBox.Critical)
            error_window.setWindowTitle("File Error")
            error_window.setInformativeText(error_message)
            error_window.setStandardButtons(QtWidgets.QMessageBox.Ok)
            error_window.exec_()
        else:
            self.finish_information_dialog()


            
    def right_click_menu(self, position):    #position是pyqt自己的參數
        
        index = self.ui_dialog.tableWidget_file_list.indexAt(position)
        
        if index.isValid() and index.column() == 0:
            right_menu = QMenu()
            
            # 檢查當前是否選擇了多行
            selected_rows = self.ui_dialog.tableWidget_file_list.selectionModel().selectedRows()
            
            if len(selected_rows) > 1:
                delete_action = right_menu.addAction("Delete Selected")
            else:
                delete_action = right_menu.addAction("Delete")
   
            # 连接菜单项的信号到相应的槽函数
            delete_action.triggered.connect(self.delete_item)
            
            
            # 在指定位置显示菜单
            right_menu.exec_(self.ui_dialog.tableWidget_file_list.viewport().mapToGlobal(position))
            
        
        
            
    def delete_item(self):
        # 獲取所有選擇的行
        selected_rows = self.ui_dialog.tableWidget_file_list.selectionModel().selectedRows()
        if not selected_rows:
            return  # 如果沒有選擇行則不執行
        
        # 逆序刪除選中的行，避免行數改變引起問題
        for index in sorted(selected_rows, reverse=True):
            row = index.row()
            item_filename = self.ui_dialog.tableWidget_file_list.item(row, 0)
            item_fileextention = self.ui_dialog.tableWidget_file_list.item(row, 1)

            if item_filename:
                file_name_in_row = item_filename.text()
                file_extention_in_row = item_fileextention.text()
                
                
                file_basename = str(file_name_in_row) + "." + str(file_extention_in_row)
         
                if file_basename in self.files_input_detail:
                    del self.files_input_detail[file_basename]
                    self.ui_dialog.tableWidget_file_list.removeRow(row)
                    
        if len(self.files_input_detail) > 0:
            self.ui_dialog.pushButton_convert.setEnabled(True)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(255, 255, 255); background-color: rgb(60, 60, 60)")
        else:
            self.ui_dialog.pushButton_convert.setEnabled(False)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(100, 100, 100); background-color: rgb(30, 30, 30)")
        
        
            
    def finish_information_dialog(self):
        info_box = QtWidgets.QMessageBox()
        info_box.setIcon(QtWidgets.QMessageBox.Information)
        info_box.setWindowTitle("Convertion Completed")
        info_box.setText("Finished")
        info_box.setStandardButtons(QtWidgets.QMessageBox.Ok)
        info_box.exec_()      
    
    
    def clear_all(self):
        self.ui_dialog.tableWidget_file_list.setRowCount(0)
        self.files_input_detail = {}
        
        if len(self.files_input_detail) > 0:
            self.ui_dialog.pushButton_convert.setEnabled(True)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(255, 255, 255); background-color: rgb(60, 60, 60)")
        else:
            self.ui_dialog.pushButton_convert.setEnabled(False)
            self.ui_dialog.pushButton_convert.setStyleSheet("font: bold 14pt \"微軟正黑體\"; color:rgb(100, 100, 100); background-color: rgb(30, 30, 30)")
        
