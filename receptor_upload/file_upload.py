# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:45:12 2024

@author: Xhamrock Studio
"""

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QHeaderView, QTableWidgetItem, QDialog, QRadioButton, QButtonGroup, QWidget, QVBoxLayout, QMenu, QCheckBox, QHBoxLayout
from PyQt5.QtCore import QTimer, QProcess, Qt, QObject, QThread, pyqtSignal


import os
import time

from MergeonDock.receptor_upload import rec_prepare_detect
from MergeonDock.progress_window import ProgressWindow
from MergeonDock.error_report import ErrorWindow


class TaskWorker(QObject):
    progress_changed = pyqtSignal(int)
    task_finished_signal = pyqtSignal()
    set_label_text_signal = pyqtSignal(str)
    show_error_signal = pyqtSignal(dict)  # ✅ 新增訊號，讓主執行緒顯示錯誤訊息
    process_stdoutput_signal = pyqtSignal(str)
    process_error_stdoutput_signal = pyqtSignal(str)
    
    def __init__(self, task_function, task_args_list):
        """
        初始化 TaskWorker
        :param task_function: 處理每個任務的函數
        :param task_args_list: 每個任務的參數列表
        """
        super().__init__()
        self.task_function = task_function
        self.task_args_list = task_args_list
        
        self.current_task = 0
        
        
        self.stdoutput_content = ""
        self.error_stdoutput_content = ""
        self.full_report = {}
        
        self.process_stdoutput_signal.connect(self.stdoutput_log_collect)
        self.process_error_stdoutput_signal.connect(self.error_stdoutput_log_collect)
        
    def stdoutput_log_collect(self, std_log_txt):
        self.stdoutput_content += std_log_txt + "\n"
        
    def error_stdoutput_log_collect(self, err_log_txt):
        self.error_stdoutput_content += err_log_txt + "\n"
        
    def run(self): 
        """
        執行所有任務
        """
        try:
            for task_args in self.task_args_list:
                # 檢查中斷請求
                if QtCore.QThread.currentThread().isInterruptionRequested():
                    print("Task interrupted.")
                    return
                
                
                try:
                    # 執行任務並傳入參數
                    self.task_function(*task_args)
                except Exception as e:
                    # 捕捉到異常，避免程式崩潰
                    task_record_name = " ".join(map(str, task_args))  #元祖轉換成字串
                    self.full_report[task_record_name] = f"{self.stdoutput_content}\n{self.error_stdoutput_content}"
                
                #清空當前輸出紀錄
                self.stdoutput_content = ""
                self.error_stdoutput_content = ""
                
                #進入第二個待處理工作
                self.current_task += 1
                progress_percentage = int((self.current_task / len(self.task_args_list)) * 100)
                self.progress_changed.emit(progress_percentage)
                

            
        finally:
            if self.full_report:
                self.show_error_signal.emit(self.full_report)  # ✅ 透過訊號傳遞錯誤資訊
                
            # 所有任務完成
            self.task_finished_signal.emit()




class Receptor_upload():
    def __init__(self, ui, pymol_process, all_parameters):
        self.ui = ui
        self.pymol_process = pymol_process  # 將 pymol_process 參數保存為類的屬性
        self.all_parameters = all_parameters
        self.progress_window = None
        self.current_edit_item = None       #用來紀錄表格中名稱新舊不同
        
        
        #UI調整
        self.ui.tableWidget_show_receptor.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ui.tableWidget_show_receptor.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.ui.tableWidget_show_receptor.resizeColumnsToContents()
        header_show_receptor = self.ui.tableWidget_show_receptor.horizontalHeader()
        header_show_receptor.setSectionResizeMode(0, QHeaderView.Stretch)               # 第0列自動伸縮
        header_show_receptor.setSectionResizeMode(1, QHeaderView.ResizeToContents)      # 第1列根據內容調整
        
        self.ui.tableWidget_show_refligands.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ui.tableWidget_show_refligands.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.ui.tableWidget_show_refligands.resizeColumnsToContents()
        header_show_refligands = self.ui.tableWidget_show_refligands.horizontalHeader()
        header_show_refligands.setSectionResizeMode(0, QHeaderView.Stretch)               
        header_show_refligands.setSectionResizeMode(1, QHeaderView.ResizeToContents)      
        
        
        # 設置表頭的左鍵點擊事件（適用於 Receptor 和 Ref Ligands 的 QTableWidget）
        header_show_receptor.sectionClicked.connect(lambda index: self.header_clicked(index, self.ui.tableWidget_show_receptor))
        header_show_refligands.sectionClicked.connect(lambda index: self.header_clicked(index, self.ui.tableWidget_show_refligands))
        
        # 初始化表頭的圖標狀態
        self.receptor_header_vis_state = False  # False 表示目前顯示分子，True 表示分子隱藏
        self.refligands_header_vis_state = False
        
        #設置右鍵菜單
        self.ui.tableWidget_show_receptor.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.tableWidget_show_receptor.customContextMenuRequested.connect(lambda position, table=self.ui.tableWidget_show_receptor: self.right_click_menu(position, table))
        self.ui.tableWidget_show_refligands.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.tableWidget_show_refligands.customContextMenuRequested.connect(lambda position, table=self.ui.tableWidget_show_refligands: self.right_click_menu(position, table))
        
        self.ui.tableWidget_show_receptor.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)      #指定用戶在點擊單元格時應選擇整行
        self.ui.tableWidget_show_receptor.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)   #允許用戶通過Ctrl 或 Shift 鍵來選擇多個項目或多個行
        
        self.ui.tableWidget_show_refligands.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.ui.tableWidget_show_refligands.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        
        #名稱變換
        self.ui.tableWidget_show_receptor.itemChanged.connect(self.receptor_item_name_changed)
        self.ui.tableWidget_show_refligands.itemChanged.connect(self.refligands_item_name_changed)
        
        #按鈕
        self.ui.pushButton_uploadreceptor.clicked.connect(self.button_upload_receptor)
        
        # 連接左鍵點擊事件
        self.ui.tableWidget_show_receptor.itemClicked.connect(self.zoom_on_click)
        self.ui.tableWidget_show_refligands.itemClicked.connect(self.zoom_on_click)
        
        
        
        
        
    
    def button_upload_receptor(self):
        receptor_path = QtWidgets.QFileDialog.getOpenFileName(None, "Choose Receptor", "", "All Supported Files (*.pdb *.pdbqt);;PDB Files (*.pdb);;PDBQT Files (*.pdbqt)" )
        
        if not receptor_path[0]:  # 如果使用者取消選擇
            return
        
        self.all_parameters.input_receptor_path = os.path.normpath(receptor_path[0]) #檔案原始完整路徑
        basename = os.path.basename(self.all_parameters.input_receptor_path)         #去除資料夾路徑(檔名+副檔名)
        self.all_parameters.input_receptor_name = os.path.splitext(basename)[0]      #去除副檔名(檔名)
        input_file_extension = os.path.splitext(basename)[1]                         #副檔名
        input_file_extension = input_file_extension.lstrip('.')
        
        
        # ✅ **如果已經有 Receptor，則先刪除舊的**
        if self.all_parameters.output_prepared_receptor_path:
            self.delete_existing_receptor()
        
        self.ui.pushButton_setgridbox.setEnabled(True)
        self.ui.pushButton_setparameter.setEnabled(True)
        self.ui.pushButton_dockingbutton.setEnabled(True)
        
        
        
        
        if os.path.isfile(self.all_parameters.input_receptor_path):
            if input_file_extension == "pdbqt":
                self.all_parameters.output_prepared_receptor_path = self.all_parameters.input_receptor_path
                self.ui.tableWidget_show_receptor.setRowCount(1)
                self.show_uploaded_receptor()
            
            elif input_file_extension != "pdbqt":
                with open(self.all_parameters.input_receptor_path, "r") as rec_file:
                    content = rec_file.readlines() 
                    
                self.HET_seq = False
                for lines in content:
                    if lines.startswith("HET   "):
                        self.HET_seq = True
                        self.open_detect_window()
                        break
                
                if self.HET_seq == False:
                    self.prepare_receptor()                               
        
        
     
        
    def prepare_receptor(self):

        # 初始化進度視窗
        self.progress_window = ProgressWindow()
        self.progress_window.show()

        
       
    
        # 建立 TaskWorker
        task_args_list = [("receptor", self.all_parameters.input_receptor_path)]
        self.worker_thread = QThread()
        self.task_worker = TaskWorker(self.run_external_process, task_args_list)
    
        # 將工作器和執行緒傳遞給 ProgressWindow
        self.progress_window.set_worker(self.worker_thread, self.task_worker)
        self.task_worker.moveToThread(self.worker_thread)
    
        # 訊號連接
        self.task_worker.progress_changed.connect(self.progress_window.set_progress_value)
        self.task_worker.set_label_text_signal.connect(self.progress_window.set_label_text)
        self.task_worker.task_finished_signal.connect(self.show_uploaded_receptor)  # ✅ 確保 UI 更新
        self.task_worker.task_finished_signal.connect(self.progress_window.process_finished)
        self.task_worker.show_error_signal.connect(self.show_error_message)
    
        # 清理執行緒
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.finished.connect(self.task_worker.deleteLater)
        self.worker_thread.started.connect(self.task_worker.run)
    
        # 啟動執行緒
        self.worker_thread.start()

     
        
        
    
    def run_external_process(self, task_type, input_path):
        """
        透過 QProcess 執行外部 AutoDock 指令
        """
        if task_type != "receptor":
            print(f"Unsupported task type: {task_type}")
            return False
    
        self.task_worker.set_label_text_signal.emit(f"Processing {os.path.basename(input_path)}...")
        
        
        # 設定轉換後的輸出檔案名稱
        output_prepared_receptor_path = os.path.normpath(
            os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.input_receptor_name}_prepared.pdbqt")
        )
    
        # 構建 AutoDock 指令
        if self.all_parameters.receptor_prepare_method == "ad4":
            if self.all_parameters.receptor_prepare_opt_switch == False:
                prepare_receptor_command = (
                    f'{self.all_parameters.autodock4_run_prepare_receptor} '
                    f'-r "{self.all_parameters.input_receptor_path}" '
                    f'-o "{output_prepared_receptor_path}"'
                )
            else:
                prepare_receptor_command = (
                    f'{self.all_parameters.autodock4_run_prepare_receptor} '
                    f'-r "{self.all_parameters.input_receptor_path}" '
                    f'-o "{output_prepared_receptor_path}" '
                    f'{self.all_parameters.autodock_prepare_receptor_custom_command}'
                )
        elif self.all_parameters.receptor_prepare_method == "meeko":
            print("Meeko support coming soon...")
            return
        
        
        
        
        
        
        # **使用 QProcess 執行指令**
        process = QtCore.QProcess()
        process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
    
        # ✅ 創建 QEventLoop 來等待結果
        event_loop = QtCore.QEventLoop()
    
        # ✅ 設定超時機制
        timeout_timer = QtCore.QTimer()
        timeout_timer.setSingleShot(True)
        timeout_timer.timeout.connect(lambda: self.on_process_timeout(process, event_loop, input_path))
    
        # **監聽 QProcess 事件**
        process.finished.connect(lambda exitCode, exitStatus: self.on_process_finished(exitCode, exitStatus, process, event_loop, input_path))
        process.errorOccurred.connect(lambda error: self.on_process_error(error, event_loop, input_path))
        process.readyReadStandardOutput.connect(lambda: self.on_process_output(process, input_path))
        process.readyReadStandardError.connect(lambda: self.on_process_output(process, input_path))
    
        # ✅ 啟動外部程式
        process.start(prepare_receptor_command)
    
        if not process.waitForStarted(5000):  # 最多等 5 秒確保啟動
            raise RuntimeError(f"⚠️ QProcess failed to start for {input_path}")
    
        # ✅ 設置超時機制
        timeout_timer.start(60000)  # **60 秒內沒結束就視為卡死**
    
        # ✅ 進入事件迴圈等待結果（但不會阻塞 UI）
        event_loop.exec_()
    
        # ✅ 檢查執行結果
        if process.exitCode() == 0:
            if os.path.exists(output_prepared_receptor_path):
                self.all_parameters.output_prepared_receptor_path = output_prepared_receptor_path
                return True
            else:
                raise RuntimeError(f"Error: Output file {output_prepared_receptor_path} not found.")
        else:
            raise RuntimeError(f"Process failed with exit code {process.exitCode()}")  # ✅ 強制拋出錯誤


    def on_process_finished(self, exitCode, exitStatus, process, event_loop, input_receptor_path):
        """當外部程式執行結束時觸發"""
        process.kill()  # 強制確保它結束
        process.waitForFinished()  # 等待確保它真的結束
        
        # 確保不管發生什麼錯誤，都結束 event_loop
        event_loop.quit()
        

    def on_process_error(self, error, event_loop, input_receptor_path):
        """當外部程式出錯時觸發"""
        
        # 確保不管發生什麼錯誤，都結束 event_loop
        event_loop.quit()
        
    
    def on_process_output(self, process, input_receptor_path):
        """即時顯示外部程式輸出"""
        output = process.readAllStandardOutput().data().decode().strip()
        error_output = process.readAllStandardError().data().decode().strip()
    
        if output:
            stdoutput_log = f"🔹 STDOUT ({input_receptor_path}): {output}"
            self.task_worker.process_error_stdoutput_signal.emit(stdoutput_log)
        if error_output:
            stderror_output_log = f"⚠️ STDERR ({input_receptor_path}): {error_output}"
            self.task_worker.process_error_stdoutput_signal.emit(stderror_output_log)
    
    
    def on_process_timeout(self, process, event_loop, input_receptor_path):
        """當外部程式超時時執行"""
        if process.state() != QtCore.QProcess.NotRunning:
            timeout_error = f"⚠️ Process timeout: {input_receptor_path} - Killing process..."
            self.task_worker.process_error_stdoutput_signal.emit(timeout_error)
            process.kill()
        event_loop.quit()  # **確保函數可以返回**
           
        
    def show_error_message(self, full_report):
        error_log_window = ErrorWindow()
        error_log_window.sorting_report_dict(full_report)
        error_log_window.exec_()
    
         
          
    
    def open_detect_window(self):
        self.open_detection = rec_prepare_detect.Receptor_sequence_detection(self.pymol_process, self.all_parameters, self)
        self.open_detection.exec_()
        
        
    def show_uploaded_receptor(self):
        if not isinstance(self.all_parameters.output_prepared_receptor_path, str):
            return  # **直接跳過函數執行，避免 TypeError**
        if os.path.exists(self.all_parameters.output_prepared_receptor_path):
            self.all_parameters.output_prepared_receptor_name = os.path.splitext(os.path.basename(self.all_parameters.output_prepared_receptor_path))[0]
            prepared_file = QTableWidgetItem(self.all_parameters.output_prepared_receptor_name)
            self.ui.tableWidget_show_receptor.setRowCount(1)
            self.ui.tableWidget_show_receptor.setItem(0, 0, prepared_file)
            self.load_file_to_pymol(self.all_parameters.output_prepared_receptor_path)
            self.pymol_process.cmd.zoom(self.all_parameters.output_prepared_receptor_name)
            
            # 創建一個 QWidget 包含 QCheckBox
            receptor_visible_widget = QWidget()
            receptor_visible_checkbox = QCheckBox()
            receptor_visible_checkbox.setChecked(True)  # 預設選中

            # 將 QCheckBox 添加到布局中
            receptor_visible_layout = QHBoxLayout()
            receptor_visible_layout.addWidget(receptor_visible_checkbox)
             
            # 調整布局，使控件居中對齊
            receptor_visible_layout.setAlignment(Qt.AlignCenter)  # 居中對齊
            receptor_visible_layout.setContentsMargins(0, 0, 0, 0)  # 設置無邊距
            
            # 將布局應用到 QWidget
            receptor_visible_widget.setLayout(receptor_visible_layout)
            self.ui.tableWidget_show_receptor.setCellWidget(0, 1, receptor_visible_widget)
 
            # 連接 QCheckBox 的信號，當狀態改變時觸發
            receptor_visible_checkbox.stateChanged.connect(lambda: self.visible_signal(receptor_visible_checkbox, self.all_parameters.output_prepared_receptor_name))
        
            
            
    
    
    def show_uploaded_ref_ligands(self):
        ref_lig_amounts = int(len(self.all_parameters.ref_prepared_ligands_path))
        if ref_lig_amounts == 0:
            return True
        self.ui.tableWidget_show_refligands.setRowCount(ref_lig_amounts)
        
        
        for ref_ligs in self.all_parameters.ref_prepared_ligands_path:
            if os.path.exists(ref_ligs):
                self.load_file_to_pymol(ref_ligs)
                ref_lig_filename = os.path.splitext(os.path.basename((ref_ligs)))[0]
                self.all_parameters.ref_prepared_ligands_name.append(ref_lig_filename)
            
        for i, ref_name in enumerate(self.all_parameters.ref_prepared_ligands_name):
            prepared_ref_lig = QTableWidgetItem(ref_name)
            self.ui.tableWidget_show_refligands.setItem(i, 0, prepared_ref_lig)
           
            
            # 創建一個 QWidget 包含 QCheckBox
            ref_ligand_visible_widget = QWidget()
            ref_ligand_visible_checkbox = QCheckBox()
            ref_ligand_visible_checkbox.setChecked(True)  # 預設選中
            
            # 將 QCheckBox 添加到布局中
            ref_ligand_visible_layout = QHBoxLayout()
            ref_ligand_visible_layout.addWidget(ref_ligand_visible_checkbox)
             
            # 調整布局，使控件居中對齊
            ref_ligand_visible_layout.setAlignment(Qt.AlignCenter)  # 居中對齊
            ref_ligand_visible_layout.setContentsMargins(0, 0, 0, 0)  # 設置無邊距
            
            # 將布局應用到 QWidget
            ref_ligand_visible_widget.setLayout(ref_ligand_visible_layout)
            self.ui.tableWidget_show_refligands.setCellWidget(i, 1, ref_ligand_visible_widget)
 
            # 連接 QCheckBox 的信號，當狀態改變時觸發
            ref_ligand_visible_checkbox.stateChanged.connect(lambda _, checkbox=ref_ligand_visible_checkbox, name=ref_name: self.visible_signal(checkbox, name))
            
        
    def delete_existing_receptor(self):
        """刪除當前的 receptor 並清理 UI"""
    
        # ✅ **從 PyMOL 中刪除舊的 receptor**
        if self.all_parameters.output_prepared_receptor_name:
            self.send_command_to_pymol(f"delete {self.all_parameters.output_prepared_receptor_name}")
    
        # ✅ **從 UI 表格中刪除**
        self.ui.tableWidget_show_receptor.setRowCount(0)
        self.ui.tableWidget_show_refligands.setRowCount(0)
    
        # ✅ **清空內部參數**
        self.all_parameters.output_prepared_receptor_name = ""
        self.all_parameters.output_prepared_receptor_path = ""
    
        # ✅ **禁用對接按鈕**
        self.ui.pushButton_dockingbutton.setEnabled(False)
        
        
        # ✅ **清理 Ref Ligands**
        if self.all_parameters.ref_prepared_ligands_name:
            for ref_ligands in self.all_parameters.ref_prepared_ligands_name:
                self.send_command_to_pymol(f"delete {ref_ligands}")
    
            self.all_parameters.ref_prepared_ligands_name = []
            self.all_parameters.ref_prepared_ligands_path = []
        
        print("✅ Existing receptor deleted.")
        
        
    
    def load_file_to_pymol(self, filepath):
        if self.pymol_process:
            try:
                self.pymol_process.cmd.load(filepath)
                
            except Exception as e:
                print("Error sending command to PyMOL:", e)
    
        
    def send_command_to_pymol(self, command):
        if self.pymol_process:
            try:
                self.pymol_process.cmd.do(command)
                print("Command sent to PyMOL:", command)
            except Exception as e:
                print("Error sending command to PyMOL:", e)
        else:
            print("PyMOL process not available.")

    

    def right_click_menu(self, position, table):    #position是pyqt自己的參數
        table = table
        index = table.indexAt(position)

        if index.isValid() and index.column() == 0:
            # 如果右鍵點擊的行未被選中，則保持當前選中的狀態，並選中當前點擊的行
            if not table.item(index.row(), 0).isSelected():
                table.selectRow(index.row())

            # 檢查當前是否選擇了多行
            selected_rows = table.selectionModel().selectedRows()

            right_menu = QMenu()
            if len(selected_rows) > 1:
                delete_action = right_menu.addAction("Delete Selected")
            else:
                delete_action = right_menu.addAction("Delete")
                
            rename_action = right_menu.addAction("Rename")
            
             
            # 连接菜单项的信号到相应的槽函数
            delete_action.triggered.connect(lambda: self.delete_item(table))
            rename_action.triggered.connect(lambda: self.rename_item(table, index.row()))
            
            # 在指定位置显示菜单
            if table == self.ui.tableWidget_show_receptor:
                right_menu.exec_(self.ui.tableWidget_show_receptor.viewport().mapToGlobal(position))
            elif table == self.ui.tableWidget_show_refligands:
                right_menu.exec_(self.ui.tableWidget_show_refligands.viewport().mapToGlobal(position))
            
            
    def delete_item(self, table):
        table = table
        selected_rows = table.selectionModel().selectedRows()
        if not selected_rows:
            return  # 如果沒有選擇行則不執行
        
        # 逆序刪除選中的行，避免行數改變引起問題
        for index in sorted(selected_rows, reverse=True):
            row = index.row()
            item = table.item(row, 0)
            
            if item:
                if table == self.ui.tableWidget_show_receptor:
                    receptor_name_raw = item.text()
                    receptor_name = receptor_name_raw.replace(' ', '_')
                    self.send_command_to_pymol(f"delete {receptor_name}")
                    self.ui.tableWidget_show_receptor.removeRow(row)
                    self.ui.tableWidget_show_refligands.setRowCount(0)
                    
                    if receptor_name == self.all_parameters.output_prepared_receptor_name:
                        self.all_parameters.output_prepared_receptor_name = ""
                        self.all_parameters.output_prepared_receptor_path = ""
                        self.ui.pushButton_dockingbutton.setEnabled(False)
                        
                        
                        for ref_ligands in self.all_parameters.ref_prepared_ligands_name:
                            self.send_command_to_pymol(f"delete {ref_ligands}")
                        
                        self.all_parameters.ref_prepared_ligands_name = []
                        self.all_parameters.ref_prepared_ligands_path = []

                        
                    
                elif table == self.ui.tableWidget_show_refligands:
                    ref_ligand_name_in_row = item.text()
                    ref_ligand_name = ref_ligand_name_in_row.replace(' ', '_')
                    self.send_command_to_pymol(f"delete {ref_ligand_name}")
                    self.ui.tableWidget_show_refligands.removeRow(row)
             
                    if ref_ligand_name in self.all_parameters.ref_prepared_ligands_name:
                        self.all_parameters.ref_prepared_ligands_name.remove(ref_ligand_name)
                        remove_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{ref_ligand_name}.pdbqt"))
                        self.all_parameters.ref_prepared_ligands_path.remove(remove_path)
                        
                       
        
        print("Current receptor:", self.all_parameters.output_prepared_receptor_path)
        print("Current ref ligand:", self.all_parameters.ref_prepared_ligands_path)
        print("Current ref ligand name:", self.all_parameters.ref_prepared_ligands_name)
       
       
         
     
    def rename_item(self, table, row):
        item = table.item(row, 0)
        
        if item:
            table.editItem(item)
               
    
    
    def receptor_item_name_changed(self, item):         #item是pyqt自己的參數(itemChanged 信号被触发时，会自动传递给槽函数)
        new_name_raw = item.text()
        new_name = new_name_raw.replace(' ', '_')  # 獲取新名稱
        #row = item.row()  # 獲取被修改的行
   
        # 獲取舊名稱
        old_name = self.all_parameters.output_prepared_receptor_name
        
        if old_name != new_name:
            self.send_command_to_pymol(f"set_name {old_name}, {new_name}")
            try:
                if self.all_parameters.output_prepared_receptor_path:
                    self.all_parameters.output_prepared_receptor_name = new_name
                    new_file_path = os.path.join(os.path.dirname(self.all_parameters.output_prepared_receptor_path), f"{new_name}.pdbqt")
                    os.rename(self.all_parameters.output_prepared_receptor_path, new_file_path)
                    self.all_parameters.output_prepared_receptor_path = new_file_path    
            except:
                self.all_parameters.output_prepared_receptor_name = old_name
                old_file_path = os.path.join(os.path.dirname(self.all_parameters.output_prepared_receptor_path), f"{old_name}.pdbqt")
                os.rename(self.all_parameters.output_prepared_receptor_path, old_file_path)
                self.all_parameters.output_prepared_receptor_path = old_file_path 
                
                error_message = "Warning message:"
                error_window = QtWidgets.QMessageBox()
                error_window.setIcon(QtWidgets.QMessageBox.Critical)
                error_window.setWindowTitle("Command Execution Error")
                error_window.setText(f"File name could not change.\nIt will still show in {new_name} but the exist file won't change, it will still be {old_name}")
                error_window.setInformativeText(error_message)
                error_window.setStandardButtons(QtWidgets.QMessageBox.Ok)
                error_window.exec_()  
           
                    
                        
            
            
    def refligands_item_name_changed(self, item):
        new_name_raw = item.text()
        new_name = new_name_raw.replace(' ', '_')  # 獲取新名稱
        row = item.row()  # 獲取被修改的行
        
        # 獲取舊名稱
        old_name = self.all_parameters.ref_prepared_ligands_name[row]
        
        # 確保新舊名稱不一樣才執行更改
        if old_name != new_name: 
          try:
              # 更新內部數據結構中的名稱
              self.all_parameters.ref_prepared_ligands_name[row] = new_name
              
              # 構造新的文件路徑
              new_file_path = os.path.join(os.path.dirname(self.all_parameters.ref_prepared_ligands_path[row]), f"{new_name}.pdbqt")
              
              if os.path.exists(new_file_path):
                  error_message = "Warning message:"
                  error_window = QtWidgets.QMessageBox()
                  error_window.setIcon(QtWidgets.QMessageBox.Critical)
                  error_window.setWindowTitle("File name exist.")
                  error_window.setText("File name was already exist. Please change another name.")
                  error_window.setInformativeText(error_message)
                  error_window.setStandardButtons(QtWidgets.QMessageBox.Ok)
                  error_window.exec_()
              else:
                  self.send_command_to_pymol(f"set_name {old_name}, {new_name}")
                  # 重命名文件
                  os.rename(self.all_parameters.ref_prepared_ligands_path[row], new_file_path)
                  
                  # 更新內部數據結構中的文件路徑
                  self.all_parameters.ref_prepared_ligands_path[row] = new_file_path
              
             
          
          except Exception as e:
              print("Error during renaming:", e)
              
              # 如果重命名失敗，還原回舊名稱
              self.all_parameters.ref_prepared_ligands_name[row] = old_name
              old_file_path = os.path.join(os.path.dirname(self.all_parameters.ref_prepared_ligands_path[row]), f"{old_name}.pdbqt")
              os.rename(self.all_parameters.ref_prepared_ligands_path[row], old_file_path)
              self.all_parameters.ref_prepared_ligands_path[row] = old_file_path 
              
              # 顯示錯誤訊息
              error_message = "Warning message:"
              error_window = QtWidgets.QMessageBox()
              error_window.setIcon(QtWidgets.QMessageBox.Critical)
              error_window.setWindowTitle("Command Execution Error")
              error_window.setText(f"File name could not change.\nIt will still show in {new_name} but the existing file won't change, it will still be {old_name}")
              error_window.setInformativeText(error_message)
              error_window.setStandardButtons(QtWidgets.QMessageBox.Ok)
              error_window.exec_()
            
        
    
    def visible_signal(self, checkbox, name):
        checkbox = checkbox
        name_raw = name
        if " " in name_raw:
            name = name_raw.replace(" ", "_")
        else:
            name = name_raw

        # 根據 QCheckBox 的狀態來發送不同的信號
        if checkbox.isChecked():
            self.pymol_process.cmd.enable(name)
            # 發送開啟的信號（ON）
        else:
            self.pymol_process.cmd.disable(name)
            # 發送關閉的信號（OFF）
    
    
    def zoom_on_click(self, item):
        # 根據點擊的內容執行 PyMOL 的 zoom 指令
        object_name = item.text()
    
        # 檢查是否需要將空格替換為底線
        if " " in object_name:
            object_name = object_name.replace(" ", "_")
    
        # 執行 PyMOL 的 zoom 命令
        self.pymol_process.cmd.zoom(object_name)
        
    
    def header_clicked(self, index, table):
        if index == 1:
            if table == self.ui.tableWidget_show_receptor:
                self.receptor_header_vis_state = not self.receptor_header_vis_state
                new_state = self.receptor_header_vis_state
                # 更改表頭的圖標
                self.update_header_icon(table, new_state)

                # 遍歷所有行的 Checkbox，並設置其狀態
                for row in range(table.rowCount()):
                    checkbox_widget = table.cellWidget(row, 1)  # 假設 Checkbox 在第 1 列
                    if checkbox_widget:
                        checkbox = checkbox_widget.findChild(QCheckBox)
                        if checkbox:
                            checkbox.setChecked(not new_state)  # 切換狀態 (開關顯示)

            elif table == self.ui.tableWidget_show_refligands:
                self.refligands_header_vis_state = not self.refligands_header_vis_state
                new_state = self.refligands_header_vis_state
                # 更改表頭的圖標
                self.update_header_icon(table, new_state)

                # 遍歷所有行的 Checkbox，並設置其狀態
                for row in range(table.rowCount()):
                    checkbox_widget = table.cellWidget(row, 1)  # 假設 Checkbox 在第 1 列
                    if checkbox_widget:
                        checkbox = checkbox_widget.findChild(QCheckBox)
                        if checkbox:
                            checkbox.setChecked(not new_state)  # 切換狀態 (開關顯示)

    def update_header_icon(self, table, state):
        header_item = QTableWidgetItem()
        
        if state == False:
            header_item.setText("👁️")  # 當顯示時，設置表頭為「眼睛」圖案
        elif state == True:
            header_item.setText("︶")  # 當隱藏時，設置表頭為「隱藏」圖案

        # 假設你要更新的是第二列（第1列，因為索引從0開始）
        table.setHorizontalHeaderItem(1, header_item)
        
        
         
        
            
 

        
        
        
        
        
