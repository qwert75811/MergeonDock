# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:16:46 2025

@author: Xhamrock Studio
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog

from MergeonDock import error_report_ui




class ErrorWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.error_ui = error_report_ui.Ui_error_Dialog()
        self.error_ui.setupUi(self)
        
        # 連接按鈕事件
        self.error_ui.error_close_pushButton.clicked.connect(self.close_window)
        self.error_ui.error_path_previous_pushButton.clicked.connect(self.show_previous_error)
        self.error_ui.error_path_next_pushButton.clicked.connect(self.show_next_error)
        
        # 儲存錯誤報告字典
        self.report_dict = {}
        self.task_names = []  # 存放 task_name 的列表
        self.current_index = 0  # 當前選擇的 index
        
        # 監聽下拉選單變化
        self.error_ui.error_path_comboBox.currentIndexChanged.connect(self.update_error_display)
        
        
    def sorting_report_dict(self, report_dict):
        self.report_dict = report_dict
        self.task_names = list(report_dict.keys())  # 取得所有 task 名稱
        
        # 清空下拉選單並重新加入 task_name
        
        self.error_ui.error_path_comboBox.clear()
        self.error_ui.error_path_comboBox.addItems(self.task_names)
        
        # 預設顯示第一個錯誤內容
        if self.task_names:
            self.current_index = 0
            self.error_ui.error_path_comboBox.setCurrentIndex(self.current_index)
            
            
    def update_error_display(self):
        """
        根據當前選擇的 task_name 顯示對應錯誤資訊
        """
        selected_task = self.error_ui.error_path_comboBox.currentText()
        if selected_task in self.report_dict:
            self.error_ui.error_log_textBrowser.setPlainText(self.report_dict[selected_task])
        else:
            self.error_ui.error_log_textBrowser.setPlainText("Error log not found.")

    def show_previous_error(self):
        """切換到前一個錯誤"""
        if self.current_index > 0:
            self.current_index -= 1
            self.error_ui.error_path_comboBox.setCurrentIndex(self.current_index)

    def show_next_error(self):
        """切換到下一個錯誤"""
        if self.current_index < len(self.task_names) - 1:
            self.current_index += 1
            self.error_ui.error_path_comboBox.setCurrentIndex(self.current_index)
        
    
    def close_window(self):
        self.close()