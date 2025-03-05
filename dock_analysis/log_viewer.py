# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:23:27 2024

@author: Xhamrock Studio
"""



from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog, QWidget, QGridLayout, QTextBrowser, QPushButton

import os

from MergeonDock.dock_analysis import log_viewer_ui


class Log_viewer(QtWidgets.QDialog):
    def __init__(self, dir_name, log_path):
        super().__init__()
        self.log_viewer_ui = log_viewer_ui.Ui_log_viewer_Form()
        self.log_viewer_ui.setupUi(self)
 
        self.ui_setting()
        
        # 記錄是否使用了初始的默認 tab
        self.initial_tab_used = False
        
        # 設置第一個 log 到默認 tab 中
        self.load_initial_log(dir_name, log_path)
        
            
    def ui_setting(self):
        # 調整 windowFlags 屬性
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowMinMaxButtonsHint | QtCore.Qt.WindowCloseButtonHint)
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowContextHelpButtonHint)
    
    def load_initial_log(self, dir_name, log_path):
        """在默認的初始 tab 中設置 log 內容"""
        if not self.initial_tab_used:
            # 設置默認 tab 的內容
            self.log_viewer_ui.tabWidget.setTabText(0, dir_name)  # 設置默認 tab 的標題

            # 讀取並顯示 log 檔案內容到初始 tab 中
            with open(log_path, "r") as file:
                log_content = file.read()
                self.log_viewer_ui.textBrowser.setText(log_content)

            # 標記默認 tab 已被使用
            self.initial_tab_used = True
        else:
            # 如果默認 tab 已經被使用過，創建新的 tab
            self.add_log(dir_name, log_path)
        
        self.log_viewer_ui.pushButton_close.clicked.connect(lambda: self.remove_tab(self.log_viewer_ui.tab))
    
    def add_log(self, dir_name, log_path):
        """新增一個新的 tab，並顯示 log 內容"""
        # 創建新 `tab` 的 `QWidget`
        new_tab_widget = QWidget()
        
        # 創建新的佈局
        layout = QGridLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setVerticalSpacing(3)
        
        # 創建新的 QTextBrowser 和關閉按鈕
        new_text_browser = QTextBrowser()
        new_text_browser.setLayoutDirection(QtCore.Qt.LeftToRight)
        new_text_browser.setStyleSheet("background-color: rgb(255, 255, 255);")
        

        new_close_button = QPushButton("Close")
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(new_close_button.sizePolicy().hasHeightForWidth())
        new_close_button.setSizePolicy(sizePolicy)
        new_close_button.setLayoutDirection(QtCore.Qt.RightToLeft)
        new_close_button.setStyleSheet("background-color: rgb(255, 255, 255);")

        
        # 設定 log 內容到 textBrowser 中
        with open(log_path, "r") as file:
            log_content = file.read()
            new_text_browser.setText(log_content)

        # 添加部件到佈局中
        layout.addWidget(new_text_browser)
        layout.addWidget(new_close_button)
        
        # 為新 tab 設定佈局
        new_tab_widget.setLayout(layout)

        # 添加新 tab 到 tabWidget 中，並設定標籤名稱為 dir_name
        self.log_viewer_ui.tabWidget.addTab(new_tab_widget, dir_name)
        self.log_viewer_ui.tabWidget.setCurrentWidget(new_tab_widget)

        # 連接新按鈕的信號，移除該 tab
        new_close_button.clicked.connect(lambda: self.remove_tab(new_tab_widget))

    
        
    
    def remove_tab(self, widget):
        """移除特定的 tab，當只剩一個時關閉整個視窗"""
        tab_count = self.log_viewer_ui.tabWidget.count()  # 取得當前 tab 數量
        index = self.log_viewer_ui.tabWidget.indexOf(widget)  # 查找 widget 對應的 tab 索引位置
    
        # 如果找到對應的 tab
        if index != -1:
            if tab_count > 1:
                # 如果 tab 數量大於 1，則只移除該 tab
                self.log_viewer_ui.tabWidget.removeTab(index)
            else:
                # 如果只剩一個 tab，關閉整個視窗
                self.close()
               
    
    
        