# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 09:55:53 2025

@author: Xhamrock Studio
"""


from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog
from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QDesktopServices

from MergeonDock.menu import about_ui



class Menu_help_about():
    def __init__(self, ui, all_parameters):
        self.ui = ui
        self.all_parameters = all_parameters
        self.ui.actionAbout.triggered.connect(self.open_own_window)
        
        
    def open_own_window(self):
        self.own_window = About()
        self.own_window.exec_()


class About(QDialog):
    def __init__(self):
        super().__init__()
        self.ui_widget = about_ui.Ui_Form_about()
        self.ui_widget.setupUi(self)
        
        self.ui_widget.textBrowser.setOpenExternalLinks(False)  # 防止內部處理
        self.ui_widget.textBrowser.anchorClicked.connect(self.open_url)
        
        # 儲存初始 HTML，當點擊後可恢復顯示
        self.original_html = self.ui_widget.textBrowser.toHtml()
        
        self.ui_widget.pushButton.clicked.connect(self.close_window)
        
    def open_url(self, url):
        """使用系統預設瀏覽器開啟超連結"""
        QDesktopServices.openUrl(url)
        
        # 重新設置原始 HTML，防止視窗變白
        QtCore.QTimer.singleShot(100, lambda: self.ui_widget.textBrowser.setHtml(self.original_html))
        
    
    
    def close_window(self):
        self.close()