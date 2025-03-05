# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:29:45 2024

@author: Xhamrock Studio
"""

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QHeaderView, QTableWidgetItem, QDialog, QRadioButton, QButtonGroup, QWidget, QVBoxLayout
from PyQt5 import QtCore
from PyQt5.QtCore import QTimer, QProcess, Qt, QObject, QThread, pyqtSignal

import os, re

from MergeonDock.receptor_upload import rec_prepare_detect_ui
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



class Receptor_sequence_detection(QDialog):
    def __init__(self, pymol_process, all_parameters, receptor_upload_instance):
        super().__init__()
        self.ui_detection = rec_prepare_detect_ui.Ui_rec_prepare_detect()
        self.ui_detection.setupUi(self)
        self.pymol_process = pymol_process
        self.all_parameters = all_parameters
        self.receptor_upload_instance = receptor_upload_instance
        
        #不用class而用def在同一個class的寫法如下
        #self.recptor_detection = rec_prepare_detect_ui.Ui_rec_prepare_detect()
        #self.recptor_detection_window = QDialog()
        #self.recptor_detection.setupUi(self.recptor_detection_window)
        #self.recptor_detection_window.show()
        
        if self.all_parameters.input_receptor_path is not None:
            self.load_file()
               
    def load_file(self):
        with open(self.all_parameters.input_receptor_path, "r") as rec_file:
            self.full_content = rec_file.read()
            self.full_line_content = self.full_content.splitlines()
             
        self.HET_residue_name = []
        self.HET_chain_id = []
        self.HET_residue_num = []

        #Sequence\Chain部分
        HET_pattern = re.compile(r'HET\s+([A-Z0-9]+)\s+([A-Z])\s*(\d+)\s+(\d+)')
        
        catch_HET_info = HET_pattern.findall(self.full_content)
        for HET_catch in catch_HET_info:
            residue_name = HET_catch[0]
            chain_id = HET_catch[1]
            residue_num = HET_catch[2]
            atom_amounts = HET_catch[3]
            
            self.HET_residue_name.append(residue_name)
            self.HET_chain_id.append(chain_id)
            self.HET_residue_num.append(residue_num)
            
                
        #Description部分   
        HET_desc = {}
        for key in self.HET_residue_name:
            HET_desc[key] = ""
        
        HET_desc_pattern = re.compile(r'HETNAM\s+([A-Z0-9]+)\s+(.+)')
        catch_desc_info = HET_desc_pattern.findall(self.full_content)
        
        for desc_catch in catch_desc_info:
            residue_name = desc_catch[0]
            description = desc_catch[1].strip()
            if residue_name in self.HET_residue_name:
                HET_desc[residue_name] = description
        
        
        
                        
        #放入表格
        self.ui_detection.tableWidget_sequence_detect.setRowCount(len(self.HET_residue_name))
        for index_seq in range(len(self.HET_residue_name)):
            self.ui_detection.tableWidget_sequence_detect.setItem(index_seq, 0, QTableWidgetItem(self.HET_residue_name[index_seq]))  
            self.ui_detection.tableWidget_sequence_detect.setItem(index_seq, 1, QTableWidgetItem(f"{self.HET_chain_id[index_seq]} | {self.HET_residue_num[index_seq]}"))
            if HET_desc[self.HET_residue_name[index_seq]]:
                self.ui_detection.tableWidget_sequence_detect.setItem(index_seq, 2, QTableWidgetItem(HET_desc[self.HET_residue_name[index_seq]])) 
        
        
        self.radio_button_group = [] # 用于存储每行的按钮组
        for row in range(len(self.HET_residue_name)):
            # 为每一行创建一个新的按钮组
            self.buttongroup = QButtonGroup(self)
            self.radio_button_group.append(self.buttongroup)
            self.remove_radiobutton = None
                        
            for col in range(3, 6):
                temp_radiobutton = QRadioButton()
                self.buttongroup.addButton(temp_radiobutton)  # 将单选按钮添加到对应行的按钮组中，实现互斥
                temp_widget = QWidget()
                temp_layout = QVBoxLayout()
                temp_layout.addWidget(temp_radiobutton)
                temp_layout.setAlignment(temp_radiobutton, Qt.AlignCenter)
                temp_layout.setContentsMargins(0, 0, 0, 0)
                temp_widget.setLayout(temp_layout)
                
                self.ui_detection.tableWidget_sequence_detect.setCellWidget(row, col, temp_widget)
                
                remove_radiobutton = temp_radiobutton
            
            if remove_radiobutton:
                remove_radiobutton.setChecked(True)
        
        
        self.ui_detection.pushButton_Abort.clicked.connect(self.abort_button)
        self.ui_detection.pushButton_Skip_preparation.clicked.connect(self.skip_preparation)
        self.ui_detection.pushButton_Contiune.clicked.connect(self.continue_)
        
        
        
    def abort_button(self):
        self.close()
    
    def skip_preparation(self):
        self.receptor_upload_instance.show_uploaded_receptor()
        self.close()    
        
    def continue_(self):
        #初始化
        self.ref_ligands_lists_for_preparation = []
        self.record_delete = []
        self.ref_ligands_detail = {}
        self.all_parameters.ref_prepared_ligands_name = []
        self.all_parameters.ref_prepared_ligands_path = []
        
        
        #進度條
        self.prepare_progress_window = ProgressWindow()
        self.prepare_progress_window.show()
        
        # 建立 task_args_list
        task_args_list = []
        
        
        # 建立 TaskWorker
        self.worker_thread = QThread()
        self.task_worker = TaskWorker(self.run_external_process, task_args_list)
        
        # 將工作器和執行緒傳遞給 ProgressWindow
        self.prepare_progress_window.set_worker(self.worker_thread, self.task_worker)
        
        self.task_worker.moveToThread(self.worker_thread)
        
        # 連接訊號
        self.task_worker.progress_changed.connect(self.prepare_progress_window.set_progress_value)
        self.task_worker.set_label_text_signal.connect(self.prepare_progress_window.set_label_text)
        self.task_worker.task_finished_signal.connect(self.receptor_upload_instance.show_uploaded_receptor)
        self.task_worker.task_finished_signal.connect(self.receptor_upload_instance.show_uploaded_ref_ligands)
        self.task_worker.task_finished_signal.connect(self.prepare_progress_window.process_finished)
        self.task_worker.show_error_signal.connect(self.show_error_message)

        # 將執行緒完成時的刪除動作與信號連接
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.finished.connect(self.task_worker.deleteLater)
        
        self.worker_thread.started.connect(self.task_worker.run)
        
        
        #收集資訊
        for row, group in enumerate(self.radio_button_group):
          selected_button = group.checkedButton() # 使用QButtonGroup的checkedButton方法来找出哪个按钮被选中
          if selected_button:
              # 在每个group中查找选中的按钮，并确定它是哪一列
              for col in range(3, 6):
                  button = self.ui_detection.tableWidget_sequence_detect.cellWidget(row, col).layout().itemAt(0).widget()
                  
                  if button == selected_button: # 確認是否為選中的按鈕
                 
                      if col == 3:
                          self.ref_ligands_preparation(row)
                          self.record_delete.append((self.HET_chain_id[row], self.HET_residue_num[row]))
                      elif col  == 4:
                          print(f"Preserve {self.HET_residue_name[row]}:{self.HET_chain_id[row]}|{self.HET_residue_num[row]}")
                          
                      else:
                          self.record_delete.append((self.HET_chain_id[row], self.HET_residue_num[row]))
        
 
        
        pattern = re.compile(r'HETATM\s*\d+\s+[A-Z0-9]{1,3}\s*([A-Z0-9]{1,4})\s+([A-Z])\s+(\d+)')
        
        content_to_keep = []
        for line in self.full_line_content:
            match = pattern.match(line)
            if match:
                residue_name = match.group(1) # 0 是完整符合pattern的字串, 1開始是括號內的
                chain_id = match.group(2)
                residue_seq = match.group(3)
                if (chain_id, residue_seq) not in self.record_delete:
                    content_to_keep.append(line)
            else:
                content_to_keep.append(line)
        
        modified_content = "\n".join(content_to_keep) + "\n"
     
        self.task_worker.set_label_text_signal.emit("Checking........Finished \nPacking info........Finished \nGenerating temp pdb File........Running")  
        self.task_worker.progress_changed.emit(50)
                          
        
        #輸出設定後的pdb暫存檔
        temp_output_receptor = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.input_receptor_name}_temp.pdb"))
        with open(temp_output_receptor, "w") as file:
            file.writelines(modified_content)  
        
        
        self.task_worker.set_label_text_signal.emit("Checking........Finished \nPacking info........Finished \nGenerating temp pdb File........Finished \nPreparing Receptor........Running")
        self.task_worker.progress_changed.emit(70)
                    
        
        #確定暫存pdb檔案後存入待處理參數表
        if os.path.exists(temp_output_receptor):
            task_args_list.append(("receptor", temp_output_receptor))
            

        if self.ref_ligands_lists_for_preparation != []:
            for unprepared_ref_ligands in self.ref_ligands_lists_for_preparation:
                task_args_list.append(("ref_ligand", unprepared_ref_ligands))
                
        
        

        self.worker_thread.start()
        

        self.close()


    def ref_ligands_preparation(self, selected_row):
        #實際位置座標的字典
        selected_seq = self.HET_residue_name[selected_row]
        selected_chain = self.HET_chain_id[selected_row]
        selected_num = self.HET_residue_num[selected_row]

        pattern = re.compile(r'HETATM\s*\d+\s+[A-Z0-9]{1,3}\s*([A-Z0-9]{1,4})\s+([A-Z])\s+(\d+)')
            
        content_of_reflig_text = []
        for line in self.full_line_content:
            match = pattern.match(line)
            if match:
                residue_name = match.group(1) # 0 是完整符合pattern的字串, 1開始是括號內的
                chain_id = match.group(2)
                residue_seq = match.group(3) 
                if chain_id == selected_chain and residue_seq == selected_num:
                    content_of_reflig_text.append(line)
                 
        self.ref_ligands_detail[selected_seq] = "\n".join(content_of_reflig_text) + "\n"

        
        #設定輸出檔案格式和路徑
        self.ref_lig_pdb_output_path = os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.input_receptor_name}_{selected_seq}_{selected_chain}_{selected_num}.pdb")
        self.output_ref_lig_path = os.path.normpath(self.ref_lig_pdb_output_path)
        #寫入原子資訊
        with open(self.output_ref_lig_path, "w") as file:
            file.write(self.ref_ligands_detail[selected_seq])
        
        self.ref_ligands_lists_for_preparation.append(self.output_ref_lig_path)


    def run_external_process(self, task_type, file_path):
        if task_type == "receptor":  
            self.receptor_pdb_to_pdbqt(file_path)  
        elif task_type == "ref_ligand":
            self.ref_ligand_to_pdbqt(file_path)  
        else:
            print(f"Unsupported task args: {task_type} {file_path}")
            return False
        
    def receptor_pdb_to_pdbqt(self, nonprepared_receptor_path):
        output_prepared_receptor_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.input_receptor_name}_prepared.pdbqt"))

        if self.all_parameters.receptor_prepare_method == "ad4":
            if self.all_parameters.receptor_prepare_opt_switch == False:
                #Optional參數未打開
                ad4_receptor_preparation_command = f'{self.all_parameters.autodock4_run_prepare_receptor} -r "{nonprepared_receptor_path}" -o "{output_prepared_receptor_path}"'
            elif self.all_parameters.receptor_prepare_opt_switch == True:
                #Optional參數打開
                ad4_receptor_preparation_command = f'{self.all_parameters.autodock4_run_prepare_receptor} -r "{nonprepared_receptor_path}" -o "{output_prepared_receptor_path}" {self.all_parameters.autodock_prepare_receptor_custom_command}'
        elif self.all_parameters.receptor_prepare_method == "meeko":
            print("coming soon...")
        
        self.task_worker.set_label_text_signal.emit("Generating receptor Files........")
        
        # **使用 QProcess 非阻塞方式，但讓函數等待結果**
        process = QtCore.QProcess()
        process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
    
        # ✅ 創建 QEventLoop 來等待結果
        event_loop = QtCore.QEventLoop()
        
        # ✅ 設定超時時間（例如 60 秒）
        timeout_timer = QtCore.QTimer()
        timeout_timer.setSingleShot(True)  #設定計時器只會觸發一次
        timeout_timer.timeout.connect(lambda: self.on_process_timeout(process, event_loop, self.all_parameters.input_receptor_name))
        
        # **監聽 QProcess 事件**
        process.finished.connect(lambda exitCode, exitStatus: self.on_process_finished(exitCode, exitStatus, process, event_loop, self.all_parameters.input_receptor_name, nonprepared_receptor_path))
        process.errorOccurred.connect(lambda error: self.on_process_error(error, event_loop, self.all_parameters.input_receptor_name))
        process.readyReadStandardOutput.connect(lambda: self.on_process_output(process, self.all_parameters.input_receptor_name))
        process.readyReadStandardError.connect(lambda: self.on_process_output(process, self.all_parameters.input_receptor_name))

    
        # ✅ 啟動外部程式
        process.start(ad4_receptor_preparation_command)

        if not process.waitForStarted(5000):  # 最多等 5 秒確保啟動
            raise RuntimeError(f"⚠️ QProcess failed to start for {self.all_parameters.input_receptor_name}")
    
        # ✅ 設置超時機制
        timeout_timer.start(60000)  # **60 秒內沒結束就視為卡死**
        
        # ✅ 進入事件迴圈等待結果（但不會阻塞 UI）
        event_loop.exec_()
    
        # ✅ 偵測最終執行結果
        if process.exitCode() == 0:
            if os.path.exists(output_prepared_receptor_path):
                if output_prepared_receptor_path != self.all_parameters.output_prepared_receptor_path:
                    self.all_parameters.output_prepared_receptor_path = output_prepared_receptor_path
                return True
            else:
                raise RuntimeError(f"Error: Output file {output_prepared_receptor_path} not found.")
        else:
            raise RuntimeError(f"Process failed with exit code {process.exitCode()}")  # ✅ 強制拋出錯誤
        

    def ref_ligand_to_pdbqt(self, nonprepared_ref_lignads_path):
        ref_ligand_filename = os.path.basename(nonprepared_ref_lignads_path)
        ref_ligand_name = os.path.splitext(ref_ligand_filename)[0]
        output_ref_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{ref_ligand_name}_ref_lig.pdbqt"))
        
        if self.all_parameters.receptor_prepare_method == "ad4":
            if self.all_parameters.ligands_prepare_opt_switch == False:
                #Optional參數未打開
                ad4_prepare_ref_ligand_command = f'{self.all_parameters.autodock4_run_prepare_ligands} -l "{nonprepared_ref_lignads_path}" -o "{output_ref_path}"'    
            elif self.all_parameters.ligands_prepare_opt_switch == True:
                #Optional參數打開
                ad4_prepare_ref_ligand_command = f'{self.all_parameters.autodock4_run_prepare_ligands} -l "{nonprepared_ref_lignads_path}" -o "{output_ref_path}" {self.all_parameters.autodock_prepare_ligands_custom_command}' 
        elif self.all_parameters.receptor_prepare_method == "meeko":
            print("coming soon...")
        
        self.task_worker.set_label_text_signal.emit("Generating ref ligands Files........")
        
        # **使用 QProcess 非阻塞方式，但讓函數等待結果**
        process = QtCore.QProcess()
        process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
    
        # ✅ 創建 QEventLoop 來等待結果
        event_loop = QtCore.QEventLoop()
        
        # ✅ 設定超時時間（例如 60 秒）
        timeout_timer = QtCore.QTimer()
        timeout_timer.setSingleShot(True)  #設定計時器只會觸發一次
        timeout_timer.timeout.connect(lambda: self.on_process_timeout(process, event_loop, f"{ref_ligand_name}_ref_lig"))
        
        # **監聽 QProcess 事件**
        process.finished.connect(lambda exitCode, exitStatus: self.on_process_finished(exitCode, exitStatus, process, event_loop, f"{ref_ligand_name}_ref_lig", nonprepared_ref_lignads_path))
        process.errorOccurred.connect(lambda error: self.on_process_error(error, event_loop, f"{ref_ligand_name}_ref_lig"))
        process.readyReadStandardOutput.connect(lambda: self.on_process_output(process, f"{ref_ligand_name}_ref_lig"))
        process.readyReadStandardError.connect(lambda: self.on_process_output(process, f"{ref_ligand_name}_ref_lig"))

    
        # ✅ 啟動外部程式
        process.start(ad4_prepare_ref_ligand_command)
        
        if not process.waitForStarted(5000):  # 最多等 5 秒確保啟動
            raise RuntimeError(f"⚠️ QProcess failed to start for {ref_ligand_name}_ref_lig")
    
        # ✅ 設置超時機制
        timeout_timer.start(60000)  # **60 秒內沒結束就視為卡死**
        
        # ✅ 進入事件迴圈等待結果（但不會阻塞 UI）
        event_loop.exec_()
    
        # ✅ 偵測最終執行結果
        if process.exitCode() == 0:
            if os.path.exists(output_ref_path):
                if output_ref_path not in self.all_parameters.ref_prepared_ligands_path:
                    self.all_parameters.ref_prepared_ligands_path.append(output_ref_path)    
                return True
            else:
                raise RuntimeError(f"Error: Output file {output_ref_path} not found.")
        else:
            raise RuntimeError(f"Process failed with exit code {process.exitCode()}")  # ✅ 強制拋出錯誤
    
    
    
    def on_process_finished(self, exitCode, exitStatus, process, event_loop, receptor_name, temp_file=None):
        """當外部程式執行結束時觸發"""
        process.kill()  # 強制確保它結束
        process.waitForFinished()  # 等待確保它真的結束
        
        # 確保不管發生什麼錯誤，都結束 event_loop
        event_loop.quit()
        
        if temp_file and os.path.exists(temp_file):
            os.remove(temp_file)

    def on_process_error(self, error, event_loop, receptor_name):
        """當外部程式出錯時觸發"""
        
        # 確保不管發生什麼錯誤，都結束 event_loop
        event_loop.quit()
        
    
    def on_process_output(self, process, receptor_name):
        """即時顯示外部程式輸出"""
        output = process.readAllStandardOutput().data().decode().strip()
        error_output = process.readAllStandardError().data().decode().strip()
    
        if output:
            stdoutput_log = f"🔹 STDOUT ({receptor_name}): {output}"
            self.task_worker.process_error_stdoutput_signal.emit(stdoutput_log)
        if error_output:
            stderror_output_log = f"⚠️ STDERR ({receptor_name}): {error_output}"
            self.task_worker.process_error_stdoutput_signal.emit(stderror_output_log)
    
    
    def on_process_timeout(self, process, event_loop, receptor_name):
        """當外部程式超時時執行"""
        if process.state() != QtCore.QProcess.NotRunning:
            timeout_error = f"⚠️ Process timeout: {receptor_name} - Killing process..."
            self.task_worker.process_error_stdoutput_signal.emit(timeout_error)
            process.kill()
        event_loop.quit()  # **確保函數可以返回**
           
        
    def show_error_message(self, full_report):
        error_log_window = ErrorWindow()
        error_log_window.sorting_report_dict(full_report)
        error_log_window.exec_()
    
    
    
        
    
       
        
            
     
    
    
  
                    
    
     
        
    
        
        
        
        
        
                        
        
                        
                        
