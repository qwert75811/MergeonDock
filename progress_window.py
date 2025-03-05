# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:48:28 2024

@author: Xhamrock Studio
"""


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog

from MergeonDock import progress_window_ui


          
            

class ProgressWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.progress_ui = progress_window_ui.Ui_Dialog_progress_window()
        self.progress_ui.setupUi(self)
        
        
        self.progress_ui.pushButton_cancel.clicked.connect(self.cancel_process)

        self.worker_thread = None  # 用於追踪工作執行緒
        self.worker = None         # 用於追踪工作器
        
    def set_worker(self, worker_thread, worker):
        """
        設置工作執行緒和工作器
        """
        self.worker_thread = worker_thread
        self.worker = worker
        
    def set_label_text(self, text):
        self.progress_ui.label_state_progress.setText(text)
    
    def set_progress_value(self, value):
        self.progress_ui.progressBar.setValue(value)
    
    def process_finished(self):
        """
        結束進度條並清理資源
        """
        self.set_label_text("Completed")
        self.set_progress_value(100)
    
        self.cleanup()
        self.close()
        

    def cancel_process(self):
        """
        用於取消工作並清理資源
        """
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.requestInterruption()  # 嘗試中斷執行緒
            self.worker_thread.quit()
            self.worker_thread.wait()
            
        self.cleanup()
        self.close()
    
    def cleanup(self):
        """
        清理執行緒和工作器資源
        """
        if self.worker and not self.worker_thread.isFinished():
            self.worker.deleteLater()
            self.worker = None  # 防止多次調用
    
        if self.worker_thread:
            self.worker_thread.quit()
            self.worker_thread.wait()  # 等待執行緒完全結束
            self.worker_thread.deleteLater()
            self.worker_thread = None  # 防止多次調用

    
    





        
        
        




