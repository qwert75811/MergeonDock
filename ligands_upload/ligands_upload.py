# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:35:24 2024

@author: Xhamrock Studio
"""

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QHeaderView, QTableWidgetItem, QMenu, QWidget, QCheckBox, QHBoxLayout
from PyQt5.QtCore import QTimer, QProcess, Qt, QObject, QThread, pyqtSignal
from PyQt5 import QtCore
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

import os
import time
import re

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



class Ligands_upload():
    def __init__(self, ui, pymol_process, all_parameters):
        self.ui = ui
        self.pymol_process = pymol_process  # 將 pymol_process 參數保存為類的屬性
        self.all_parameters = all_parameters
         
        
        #UI調整
        self.ui.tableWidget_show_ligands.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        header_show_ligands = self.ui.tableWidget_show_ligands.horizontalHeader()
        header_show_ligands.setSectionResizeMode(0, QHeaderView.Stretch)               # 第0列自動伸縮
        header_show_ligands.setSectionResizeMode(1, QHeaderView.ResizeToContents)      # 第1列根據內容調整
        
        # 設置表頭的左鍵點擊事件（適用於 Receptor 和 Ref Ligands 的 QTableWidget）
        header_show_ligands.sectionClicked.connect(lambda index: self.header_clicked(index, self.ui.tableWidget_show_ligands))

        # 初始化表頭的圖標狀態
        self.ligands_header_vis_state = False  # False 表示目前顯示分子，True 表示分子隱藏
        

        self.ui.tableWidget_show_ligands.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)      #指定用戶在點擊單元格時應選擇整行
        self.ui.tableWidget_show_ligands.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

        
        #設置右鍵菜單
        self.ui.tableWidget_show_ligands.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.tableWidget_show_ligands.customContextMenuRequested.connect(self.right_click_menu)
        
        #名稱變換
        self.ui.tableWidget_show_ligands.itemChanged.connect(self.item_name_changed)
        
        #按鈕
        self.ui.pushButton_uploadligands.clicked.connect(self.button_upload_ligands)
        
        # 連接左鍵點擊事件
        self.ui.tableWidget_show_ligands.itemClicked.connect(self.zoom_on_click)
        
    
    def button_upload_ligands(self):  
        ligands_lists_raw = QtWidgets.QFileDialog.getOpenFileNames(None, "Choose ligands", "", "All Supported Files (*.pdb *.pdbqt *.sdf);;PDB Files (*.pdb);;PDBQT Files (*.pdbqt);;SDF Files (*.sdf)")
        
        if ligands_lists_raw:
            self.prepared_ligands_dic = os.path.normpath(os.path.join(self.all_parameters.work_directory, "prepared_lignads"))       
            os.makedirs(self.prepared_ligands_dic, exist_ok=True)
            
            ligands_path_list = []
            for files_path in ligands_lists_raw[0]:
                ligands_path_list.append(files_path)
            
            self.prepare_ligands(ligands_path_list)
        else:
            print("File not found")
        
 

            
    def prepare_ligands(self, ligands_path_list):
        total_tasks = len(ligands_path_list)
        if total_tasks == 0:
            print("No ligands to prepare.")
            return
        
        # 打包參數列表
        task_args_list = [(path,) for path in ligands_path_list]
        
        # 初始化進度條窗口
        self.prepare_progress_window = ProgressWindow()
        self.prepare_progress_window.show()
    
        # 初始化 TaskWorker
        self.worker_thread = QThread()
        self.task_worker = TaskWorker(self.run_single_ligand_task, task_args_list)
        
        # 將工作器和執行緒傳遞給 ProgressWindow
        self.prepare_progress_window.set_worker(self.worker_thread, self.task_worker)
    
        # 將 TaskWorker 信號與進度條連接
        self.task_worker.progress_changed.connect(self.prepare_progress_window.set_progress_value)
        self.task_worker.set_label_text_signal.connect(self.prepare_progress_window.set_label_text)
        self.task_worker.task_finished_signal.connect(self.load_ligands_to_ui)
        self.task_worker.task_finished_signal.connect(self.prepare_progress_window.process_finished)
    
        # ✅ 讓 `show_error_signal` 連接到主執行緒的錯誤顯示函數
        self.task_worker.show_error_signal.connect(self.show_error_message)
        
        
        # 將執行緒完成時的刪除動作與信號連接
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.finished.connect(self.task_worker.deleteLater)
    
        # 啟動 TaskWorker
        self.task_worker.moveToThread(self.worker_thread)
        self.worker_thread.started.connect(self.task_worker.run)
        
        # 啟動執行緒
        self.worker_thread.start()

    
    def run_single_ligand_task(self, *args):
        if len(args) == 1:  # 單一檔案
            file_path = args[0]
            extension = os.path.splitext(file_path)[1].lstrip('.')
            if extension == "pdb":
                return self.pdb_to_pdbqt(file_path)
            elif extension == "sdf":
                return self.import_convert(extension, file_path)
        elif len(args) == 2:  # 分子資訊 (來自 SDF)
            ligand_name, molecule_info = args
            return self.convert_into_pdb(ligand_name, molecule_info)
        else:
            print(f"Unsupported task args: {args}")
            return False

         
            
    def import_convert(self, extension, file_path):
        if extension == "sdf":
            ligands_info = self.extract_sdf(file_path)  # 取得 SDF 裡面的分子資訊
            if ligands_info:
                # 將每個分子新增為子任務
                for ligand_name, molecule_info in ligands_info.items():
                    self.task_worker.task_args_list.append((ligand_name, molecule_info))
            else:
                print("Error: No ligands found in SDF.")
                return False
        elif extension in ("csv", "tsv"):
            sdf_file = self.convert_into_sdf(file_path)
            self.import_convert("sdf", sdf_file)  # 對 SDF 文件進行處理
        else:
            print(f"Unsupported file type: {extension}")
        return True
        
        
    
    def pdb_to_pdbqt(self, file_path):
        ligand_basename = os.path.basename(file_path) 
        ligand_name = os.path.splitext(ligand_basename)[0]
        
        self.output_prepared_ligands_path = os.path.normpath(os.path.join(self.prepared_ligands_dic, f"{ligand_name}.pdbqt"))
        if self.all_parameters.ligands_prepare_opt_switch == False:
            ad4_prepare_ligands = f'{self.all_parameters.autodock4_run_prepare_ligands} -l "{file_path}" -o "{self.output_prepared_ligands_path}"'
        elif self.all_parameters.ligands_prepare_opt_switch == True:
            ad4_prepare_ligands = f'{self.all_parameters.autodock4_run_prepare_ligands} -l "{file_path}" -o "{self.output_prepared_ligands_path}" {self.all_parameters.autodock_prepare_ligands_custom_command}' 


 
        # **使用 QProcess 非阻塞方式，但讓函數等待結果**
        process = QtCore.QProcess()
        process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
    
        # ✅ 創建 QEventLoop 來等待結果
        event_loop = QtCore.QEventLoop()
        
        # ✅ 設定超時時間（例如 60 秒）
        timeout_timer = QtCore.QTimer()
        timeout_timer.setSingleShot(True)  #設定計時器只會觸發一次
        timeout_timer.timeout.connect(lambda: self.on_process_timeout(process, event_loop, ligand_name))
        
        # **監聽 QProcess 事件**
        process.finished.connect(lambda exitCode, exitStatus: self.on_process_finished(exitCode, exitStatus, process, event_loop, ligand_name))
        process.errorOccurred.connect(lambda error: self.on_process_error(error, event_loop, ligand_name))
        process.readyReadStandardOutput.connect(lambda: self.on_process_output(process, ligand_name))
        process.readyReadStandardError.connect(lambda: self.on_process_output(process, ligand_name))

    
        # ✅ 啟動外部程式
        process.start(ad4_prepare_ligands)
        
        if not process.waitForStarted(5000):  # 最多等 5 秒確保啟動
            raise RuntimeError(f"⚠️ QProcess failed to start for {ligand_name}")
    
        # ✅ 設置超時機制
        timeout_timer.start(60000)  # **60 秒內沒結束就視為卡死**
        
        # ✅ 進入事件迴圈等待結果（但不會阻塞 UI）
        event_loop.exec_()
    
        # ✅ 偵測最終執行結果
        if process.exitCode() == 0:
            if os.path.exists(self.output_prepared_ligands_path):
                if self.output_prepared_ligands_path not in self.all_parameters.output_prepared_ligands_path:
                    self.all_parameters.output_prepared_ligands_path.append(self.output_prepared_ligands_path)
                return True
            else:
                raise RuntimeError(f"Error: Output file {self.output_prepared_ligands_path} not found.")
        else:
            raise RuntimeError(f"Process failed with exit code {process.exitCode()}")  # ✅ 強制拋出錯誤
        
        
    def on_process_finished(self, exitCode, exitStatus, process, event_loop, ligand_name):
        """當外部程式執行結束時觸發"""
        process.kill()  # 強制確保它結束
        process.waitForFinished()  # 等待確保它真的結束
        
        # 確保不管發生什麼錯誤，都結束 event_loop
        event_loop.quit()
        

    def on_process_error(self, error, event_loop, ligand_name):
        """當外部程式出錯時觸發"""
        
        # 確保不管發生什麼錯誤，都結束 event_loop
        event_loop.quit()
        
    
    def on_process_output(self, process, ligand_name):
        """即時顯示外部程式輸出"""
        output = process.readAllStandardOutput().data().decode().strip()
        error_output = process.readAllStandardError().data().decode().strip()
    
        if output:
            stdoutput_log = f"🔹 STDOUT ({ligand_name}): {output}"
            self.task_worker.process_error_stdoutput_signal.emit(stdoutput_log)
        if error_output:
            stderror_output_log = f"⚠️ STDERR ({ligand_name}): {error_output}"
            self.task_worker.process_error_stdoutput_signal.emit(stderror_output_log)
    
    
    def on_process_timeout(self, process, event_loop, ligand_name):
        """當外部程式超時時執行"""
        if process.state() != QtCore.QProcess.NotRunning:
            timeout_error = f"⚠️ Process timeout: {ligand_name} - Killing process..."
            self.task_worker.process_error_stdoutput_signal.emit(timeout_error)
            process.kill()
        event_loop.quit()  # **確保函數可以返回**
           
        
    def show_error_message(self, full_report):
        error_log_window = ErrorWindow()
        error_log_window.sorting_report_dict(full_report)
        error_log_window.exec_()
        
    
    
 
    def extract_sdf(self, file_path):
        sdf_content = file_path
        ligands_info = {}
        if os.path.exists(sdf_content):
            with open(sdf_content, "r") as file:
                file_content = file.read()
                
                # 分割檔案中的每個分子塊，以 '$$$$' 為界限
                molecules = re.split(r'\$\$\$\$\n', file_content)
 
                #標籤匹配規則
                tag_pattern = re.compile(r'>\s*<([^>]+)>\s*(.*?)\n', re.MULTILINE)     
                
                # 遍歷每個分子塊
                for molecule in molecules:
                    molecule_info = {"molecule_content": molecule,"tags": {}}     #小字典: 分子資訊存取
                    
                    tags = list(tag_pattern.finditer(molecule))  #針對標籤查找

                    #找到第一個標籤的值當作大字典的鍵
                    if tags:
                        first_tag = tags[0]
                        first_tag_syn = first_tag.group(1).strip() 
                        first_tag_content = first_tag.group(2).strip()
                    
                        # 匹配这个分子块中的所有标签和它们的值
                        for match in tags:
                            tag = match.group(1)
                            value = match.group(2)
                            molecule_info["tags"][tag] = value  # 把全部的标签和值存进小字典
                        
                        ligands_info[first_tag_content]= molecule_info  #把剛剛的第一個標籤的值當作大字典的鍵 把上面小字典的內容存在這
                return ligands_info
                              
        else:
            print("they are string")
            return False
            
            
    def convert_into_sdf(self, all_content):
        print("Building soon....")
        
        
    def convert_into_pdb(self, ligand_name, molecule_info):
        # 确保目录存在，如果不存在则创建
        temp_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, "temp"))       
        os.makedirs(temp_path, exist_ok=True)
        
        failed_list = {}
        
        self.task_worker.set_label_text_signal.emit(f"Converting {ligand_name}...")
        
        molecule_content = molecule_info["molecule_content"]

        # 从字符串读取SDF内容
        molecule = Chem.MolFromMolBlock(molecule_content, strictParsing=True)

        if molecule:
            fragments = self.check_and_split_fragments(molecule)
            
            # 例如保留最大的片段作為主要片段
            main_fragment = max(fragments, key=lambda frag: frag.GetNumAtoms())
            
            try:
                # 添加氢原子
                main_fragment = Chem.AddHs(main_fragment, addCoords=True)

                # 使用 MMFF 力場優化分子結構
                AllChem.EmbedMolecule(main_fragment, randomSeed=42)  # 嵌入分子的 3D 坐標
                    
                """
                    try:
                    # 嘗試 MMFF 力場優化
                        mmff_props = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant='MMFF94')
                        if mmff_props is not None:
                            optimization_result = AllChem.MMFFOptimizeMolecule(molecule, mmff_props)
                            if optimization_result != 0:
                                raise ValueError(f"MMFF optimization failed for molecule {ligand_name}")
                    except Exception as e:
                        # 如果 MMFF 力場優化失敗，嘗試使用 UFF
                        optimization_result = AllChem.UFFOptimizeMolecule(molecule)
                        if optimization_result != 0:
                            raise ValueError(f"UFF optimization failed for molecule {ligand_name}")
                """
                
                # 使用 SanitizeMol 对分子进行“消毒”
                Chem.SanitizeMol(main_fragment)

                # 计算Gasteiger电荷
                AllChem.ComputeGasteigerCharges(main_fragment)

                # 转换为PDB格式
                pdb_block = Chem.MolToPDBBlock(main_fragment)
                
                pdb_header = f"COMPND    {ligand_name}\n"
                pdb_content = pdb_header + pdb_block

                # 写入PDB文件
                output_file = os.path.normpath(os.path.join(temp_path, f"{ligand_name}.pdb"))
                with open(output_file, "w") as file:
                    file.write(pdb_content)
                
                self.task_worker.task_args_list.append((output_file,))
                
     
                
            except Exception as e:
                failed_list[ligand_name] = f"Failed to process molecule {ligand_name}: {str(e)}"
                
        else:
            failed_list[ligand_name] = f"Failed to read molecule content for {ligand_name}."
    
        return True
     
    def check_and_split_fragments(self, molecule):
        # 使用 GetMolFrags 獲取所有片段
        fragments = rdmolops.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    
        # 返回片段列表
        return fragments   



    def load_ligands_to_ui(self):
        """
        將轉換完成的 ligands 載入到 PyMOL 和 UI，並顯示進度
        """
        total_ligands = len(self.all_parameters.output_prepared_ligands_path)
        if total_ligands == 0:
            print("No ligands to load.")
            return
    
        # 預設 UI 行數
        self.ui.tableWidget_show_ligands.setRowCount(total_ligands)
    
        for i, prepared_ligand in enumerate(self.all_parameters.output_prepared_ligands_path, start=1):
            if os.path.exists(prepared_ligand):
                ligand_basename = os.path.basename(prepared_ligand)
                ligand_name = os.path.splitext(ligand_basename)[0]
                if ligand_name not in self.all_parameters.output_prepared_ligands_name:
                    self.all_parameters.output_prepared_ligands_name.append(ligand_name)
                    self.load_file_to_pymol(prepared_ligand)
    
                # 單行更新表格
                self.update_table_row(i - 1, ligand_name)
    
            # 更新進度條
            progress_percentage = int((i / total_ligands) * 100)
            self.prepare_progress_window.set_progress_value(progress_percentage)
            self.prepare_progress_window.set_label_text(f"Loading {ligand_name} ({i}/{total_ligands})")

    def update_table_row(self, row_index, lig_name):
        """
        僅更新單行數據，減少表格重繪次數
        """
        prepared_lig = QTableWidgetItem(lig_name)
        self.ui.tableWidget_show_ligands.setItem(row_index, 0, prepared_lig)
    
        # 創建 QCheckBox
        ligand_visible_widget = QWidget()
        ligand_visible_checkbox = QCheckBox()
        ligand_visible_checkbox.setChecked(True)  # 預設選中
    
        # 設置佈局
        ligand_visible_layout = QHBoxLayout()
        ligand_visible_layout.addWidget(ligand_visible_checkbox)
        ligand_visible_layout.setAlignment(Qt.AlignCenter)
        ligand_visible_layout.setContentsMargins(0, 0, 0, 0)
        ligand_visible_widget.setLayout(ligand_visible_layout)
    
        self.ui.tableWidget_show_ligands.setCellWidget(row_index, 1, ligand_visible_widget)
    
        # 信號連接
        ligand_visible_checkbox.stateChanged.connect(lambda _, checkbox=ligand_visible_checkbox, row = row_index : self.visible_signal(checkbox, self.ui.tableWidget_show_ligands.item(row, 0).text()))
        
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
            
            
            
    def right_click_menu(self, position):    #position是pyqt自己的參數
        
        index = self.ui.tableWidget_show_ligands.indexAt(position)
        
        if index.isValid() and index.column() == 0:
            right_menu = QMenu()
            
            # 檢查當前是否選擇了多行
            selected_rows = self.ui.tableWidget_show_ligands.selectionModel().selectedRows()
            
            if len(selected_rows) > 1:
                delete_action = right_menu.addAction("Delete Selected")
            else:
                delete_action = right_menu.addAction("Delete")
                
            rename_action = right_menu.addAction("Rename")
             
            
            # 连接菜单项的信号到相应的槽函数
            delete_action.triggered.connect(self.delete_item)
            rename_action.triggered.connect(lambda: self.rename_item(index.row()))
            
            # 在指定位置显示菜单
            right_menu.exec_(self.ui.tableWidget_show_ligands.viewport().mapToGlobal(position))
            
        
        
            
    def delete_item(self):
        # 獲取所有選擇的行
        selected_rows = self.ui.tableWidget_show_ligands.selectionModel().selectedRows()
        if not selected_rows:
            return  # 如果沒有選擇行則不執行
        
        # 逆序刪除選中的行，避免行數改變引起問題
        for index in sorted(selected_rows, reverse=True):
            row = index.row()
            item = self.ui.tableWidget_show_ligands.item(row, 0)

            if item:
                ligand_name_in_row = item.text()
                ligand_name = ligand_name_in_row.replace(' ', '_')
                self.send_command_to_pymol(f"delete {ligand_name}")
                self.ui.tableWidget_show_ligands.removeRow(row)
         
                if ligand_name in self.all_parameters.output_prepared_ligands_name:
                    self.all_parameters.output_prepared_ligands_name.remove(ligand_name)
                    
                    remove_path = os.path.normpath(os.path.join(self.prepared_ligands_dic, f"{ligand_name}.pdbqt"))
                    self.all_parameters.output_prepared_ligands_path.remove(remove_path)
                    
                    
                            
               
                
                    
        print("Current ligands:", self.all_parameters.output_prepared_ligands_path)
        print("Current ligands name:", self.all_parameters.output_prepared_ligands_name)
        
      
            
            
     
    def rename_item(self, row):
        item = self.ui.tableWidget_show_ligands.item(row, 0)
        if item:
            self.ui.tableWidget_show_ligands.editItem(item)
    
    
    def item_name_changed(self, item):         #item是pyqt自己的參數(itemChanged 信号被触发时，会自动传递给槽函数)
        new_name_raw = item.text()
        new_name = new_name_raw.replace(' ', '_')  # 獲取新名稱
        row = item.row()  # 獲取被修改的行
   
        # 獲取舊名稱
        old_name = self.all_parameters.output_prepared_ligands_name[row]
   
        # 確保新舊名稱不一樣才執行更改
        if old_name != new_name: 
          try:
              # 更新內部數據結構中的名稱
              self.all_parameters.output_prepared_ligands_name[row] = new_name
              
              # 構造新的文件路徑
              new_file_path = os.path.join(os.path.dirname(self.all_parameters.output_prepared_ligands_path[row]), f"{new_name}.pdbqt")
              
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
                  os.rename(self.all_parameters.output_prepared_ligands_path[row], new_file_path)
                  
                  # 更新內部數據結構中的文件路徑
                  self.all_parameters.output_prepared_ligands_path[row] = new_file_path
              
             
          
          except Exception as e:
              print("Error during renaming:", e)
              
              # 如果重命名失敗，還原回舊名稱
              self.all_parameters.output_prepared_ligands_name[row] = old_name
              old_file_path = os.path.join(os.path.dirname(self.all_parameters.output_prepared_ligands_path[row]), f"{old_name}.pdbqt")
              os.rename(self.all_parameters.output_prepared_ligands_path[row], old_file_path)
              self.all_parameters.output_prepared_ligands_path[row] = old_file_path 
              
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
            if table == self.ui.tableWidget_show_ligands:
                self.ligands_header_vis_state = not self.ligands_header_vis_state
                new_state = self.ligands_header_vis_state
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
    
    
    
# 程式邏輯流程：
# 1. 使用者點擊按鈕 (button_upload_ligands)
#    選擇上傳的分子檔案，然後呼叫 prepare_ligands(ligands_path_list) 來準備處理。

# 2. 建立 TaskWorker 並啟動 QThread
#    prepare_ligands() 會建立 TaskWorker 物件，並將其移動到 QThread 內運行。
#    TaskWorker 接受 run_single_ligand_task 作為處理函數，開始依序處理所有任務。

# 3. TaskWorker 執行任務 (run)
#    TaskWorker 會依序執行 self.task_function(*task_args)，即 run_single_ligand_task。
#    這個函數會根據檔案類型決定要執行 pdb_to_pdbqt 或 import_convert 來轉換檔案。

# 4. QProcess 啟動外部 AutoDock (pdb_to_pdbqt)
#    QProcess 非同步執行 AutoDock 轉換指令，並監聽其輸出 (on_process_output)。
#    任何錯誤或標準輸出都會透過 process_stdoutput_signal & process_error_stdoutput_signal 傳遞回 TaskWorker。

# 5. TaskWorker 接收並記錄輸出
#    process_stdoutput_signal & process_error_stdoutput_signal 會將輸出累積 (stdoutput_log_collect, error_stdoutput_log_collect)。
#    當 TaskWorker 執行失敗時，會將對應檔案的 stdout + stderr 記錄到 full_report。

# 6. 任務完成後回報結果
#    當所有任務執行完畢，TaskWorker 會發送 task_finished_signal，通知 Ligands_upload 載入完成的分子 (load_ligands_to_ui)。
#    如果有錯誤，則透過 show_error_signal 顯示錯誤訊息 (show_error_message)。

# 訊號傳遞總結：
# progress_changed → 更新進度條
# task_finished_signal → 通知 Ligands_upload 載入結果
# show_error_signal → 通知 Ligands_upload 顯示錯誤訊息
# process_stdoutput_signal / process_error_stdoutput_signal → 記錄 stdout/stderr 並存入 full_report  
            
        
      
            
        
    
    
           