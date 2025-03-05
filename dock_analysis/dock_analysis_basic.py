# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 15:56:21 2024

@author: Xhamrock Studio
"""

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QHeaderView, QTableWidgetItem, QTableWidget, QCheckBox, QWidget, QHBoxLayout, QLabel, QMessageBox, QPushButton, QComboBox, QFileDialog, QVBoxLayout, QDialog

from PyQt5.QtCore import Qt, pyqtSignal
import os
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from openbabel import openbabel
from MergeonDock.dock_analysis import log_viewer


class Analysis_results():
    def __init__(self, ui, pymol_process, all_parameters):
        self.ui = ui
        self.pymol_process = pymol_process  
        self.all_parameters = all_parameters
        
        self.ui_advance_setup()
        
        
        
        # 暫時禁用 auto_zoom
        self.pymol_process.cmd.set("auto_zoom", 0)
        
    
        
    def ui_advance_setup(self):
        self.ui.tableWidget_analysis_receptor.resizeColumnsToContents()
        self.ui.tableWidget_analysis_receptor.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        
        header_analysis_receptor = self.ui.tableWidget_analysis_receptor.horizontalHeader()
        header_analysis_receptor.setSectionResizeMode(0, QHeaderView.ResizeToContents)  # 第0列根據內容調整
        header_analysis_receptor.setSectionResizeMode(1, QHeaderView.ResizeToContents)           # 第1列根據內容調整
        header_analysis_receptor.setSectionResizeMode(2, QHeaderView.Stretch)           # 第2列自動伸縮
        
        
 
        # 連接左鍵點擊事件
        self.ui.tableWidget_analysis_receptor.itemClicked.connect(self.zoom_on_click)


        self.ui.tableWidget_analysis_ligands.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.ui.tableWidget_analysis_ligands.resizeColumnsToContents()
        header_analysis_ligands = self.ui.tableWidget_analysis_ligands.horizontalHeader()

        # 針對特定列設置不同的調整模式
        header_analysis_ligands.setSectionResizeMode(0, QHeaderView.Stretch)  
        header_analysis_ligands.setSectionResizeMode(1, QHeaderView.Stretch)           
        header_analysis_ligands.setSectionResizeMode(2, QHeaderView.ResizeToContents)         
        header_analysis_ligands.setSectionResizeMode(3, QHeaderView.ResizeToContents) 
        self.ui.tableWidget_analysis_ligands.setColumnWidth(1, 40)
        
        
        # 設置表頭的左鍵點擊事件（適用於 Receptor 和 Ref Ligands 的 QTableWidget）
        header_analysis_receptor.sectionClicked.connect(lambda index: self.header_clicked(index, self.ui.tableWidget_analysis_receptor))
        header_analysis_ligands.sectionClicked.connect(lambda index: self.header_clicked(index, self.ui.tableWidget_analysis_ligands))
        
        
        
        # 初始化表頭的圖標狀態
        self.ana_receptor_header_vis_state = False
        self.ana_ligands_header_vis_state = False  # False 表示目前顯示分子，True 表示分子隱藏
        
        self.interaction_analysis_vis_state = False
        
        
        
        header_interaction_analysis = self.ui.tableWidget_interaction_analysis.horizontalHeader()
        header_interaction_analysis.setSectionResizeMode(0, QHeaderView.ResizeToContents)  # Interaction 根據內容調整
        header_interaction_analysis.setSectionResizeMode(1, QHeaderView.Stretch)           # Atom 讓它自動延展
        header_interaction_analysis.setSectionResizeMode(2, QHeaderView.ResizeToContents)  # Residue 根據內容調整
        header_interaction_analysis.setSectionResizeMode(3, QHeaderView.ResizeToContents)  # Distance 根據內容調整
        header_interaction_analysis.setSectionResizeMode(4, QHeaderView.ResizeToContents)  # Checkbox (是否顯示)
        
        
        
        header_interaction_analysis.sectionClicked.connect(lambda index: self.header_clicked(index, self.ui.tableWidget_interaction_analysis))

        self.ui.tableWidget_interaction_analysis.setSortingEnabled(True)

        
        # 確保按鈕狀態同步
        self.update_save_function_button()
  
        
        #按鈕
        self.ui.pushButton_analysis_load.clicked.connect(self.button_load_data)
        self.ui.pushButton_analysis_search.clicked.connect(self.affinity_filter_search)
        self.ui.pushButton_analysis_reset.clicked.connect(self.affinity_filter_reset)
        
        self.ui.pushButton_analysis_save_image.clicked.connect(self.save_image_action)
        self.ui.pushButton_analysis_save_complex.clicked.connect(self.save_complex_action)
        self.ui.pushButton_analysis_save_ligand.clicked.connect(self.save_ligand_action)
        self.ui.pushButton_interaction_analysis.clicked.connect(self.toggle_interaction_analysis)
        
        self.ui.pushButton_analysis_back.clicked.connect(self.back_to_analysis_result)
        
        self.ui.pushButton_interaction_save_image.clicked.connect(self.save_image_action)
        
        self.ui.pushButton_interaction_save_complex.clicked.connect(self.save_interaction_action) 
        

        #功能暫時禁用(保留)
        self.ui.pushButton_analysis_add.setEnabled(False)
     
    
    def update_save_function_button(self):
        """
        根據表格內容來啟用或禁用 Save Image 按鈕
        """
        if self.ui.tableWidget_analysis_receptor.rowCount() > 0:
            self.ui.pushButton_analysis_save_image.setEnabled(True)  # 啟用按鈕
            self.ui.pushButton_analysis_save_complex.setEnabled(True)
            self.ui.pushButton_analysis_save_ligand.setEnabled(True)
            self.ui.pushButton_interaction_analysis.setEnabled(True)
        else:
            self.ui.pushButton_analysis_save_image.setEnabled(False)  # 禁用按鈕
            self.ui.pushButton_analysis_save_complex.setEnabled(False)
            self.ui.pushButton_analysis_save_ligand.setEnabled(False)
            self.ui.pushButton_interaction_analysis.setEnabled(False)
    
        
    def auto_load_from_dock_tab(self, cdl_path):
        cdl_path = cdl_path
        self.ui.tabWidget.setCurrentWidget(self.ui.tab_analysis)
        self.ui.stackedWidget_analysis.setCurrentWidget(self.ui.page_analysis_basic)
        self.current_directory = self.all_parameters.work_directory
        self.pymol_process.cmd.reinitialize()
        self.data_dict = {f"{cdl_path}":{"extract_data":{}, "ligand_data":{}}}
        self.ui.tableWidget_analysis_receptor.setRowCount(0)
        self.ui.tableWidget_analysis_ligands.setRowCount(0)
        self.ui.tableWidget_interaction_analysis.setRowCount(0)
        self.load_cdl_data(cdl_path)
        
        # 確保按鈕狀態同步
        self.update_save_function_button()
    
    def button_load_data(self):
        if self.ui.tableWidget_analysis_receptor.rowCount() != 0:
            question_window = QMessageBox()
            question_window.setIcon(QMessageBox.Question)
            question_window.setWindowTitle("Notice")
            question_window.setText("Do you want to remove current section for new section?")
            question_window.setStandardButtons(QMessageBox.Yes | QMessageBox.No)  # 添加 Yes 和 No 按鈕
            question_window.setDefaultButton(QMessageBox.No)  # 預設選擇 No 按鈕
            
            question_window_reply = question_window.exec_()
            
            if question_window_reply == question_window.Yes:
                data_path = QtWidgets.QFileDialog.getOpenFileName(None, "Choose cdl file", "", "cdl files (*.cdl)")
                if not data_path[0]:
                    return
                self.ui.tableWidget_analysis_receptor.setRowCount(0)
                upload_path = os.path.normpath(data_path[0])
                self.current_directory = os.path.dirname(upload_path)
                basename = os.path.basename(upload_path)
                file_extention = os.path.splitext(basename)[1]
                self.data_dict = {}
                
                if file_extention == ".cdl":
                    self.data_dict = {f"{upload_path}":{"extract_data":{}, "ligand_data":{}}}
                    self.pymol_process.cmd.reinitialize()
                    self.load_cdl_data(upload_path)
                elif file_extention == ".txt":
                    self.pymol_process.cmd.reinitialize()
                    self.load_txt_data()
                
            else:
                print("cancel")
        else:
            data_path = QtWidgets.QFileDialog.getOpenFileName(None, "Choose cdl file", "", "cdl files (*.cdl)")
            upload_path = os.path.normpath(data_path[0])
            self.current_directory = os.path.dirname(upload_path)
            basename = os.path.basename(upload_path)
            file_extention = os.path.splitext(basename)[1]
            self.data_dict = {}
            
            if file_extention == ".cdl":
                self.data_dict = {f"{upload_path}":{"extract_data":{}, "ligand_data":{}}}
                self.pymol_process.cmd.reinitialize()
                self.load_cdl_data(upload_path)
            elif file_extention == ".txt":
                self.pymol_process.cmd.reinitialize()
                self.load_txt_data()
        
        
        # 確保按鈕狀態同步
        self.update_save_function_button()
    
    
    
    def load_cdl_data(self, upload_path):
        #獲取資訊---------------------------------------------------------------------------------------------------------------
        upload_path = upload_path
        extract_info = self.data_dict[f"{upload_path}"]["extract_data"]
        result_info = self.data_dict[f"{upload_path}"]["ligand_data"]
        
        with open(upload_path, "r") as file:
            content_line = file.readlines()
            for line in content_line:
                if line.startswith("Work directory: "):
                    extract_info["extract_dirpath"] = line.split(":", 1)[1].strip()   #split(sep, maxsplit)                  
                elif line.startswith("Receptor:"):
                    extract_info["extract_receptor"] = line.split(":", 1)[1].strip()
                elif line.startswith("Ligands:"): 
                    extract_info["extract_ligands"] = eval(line.split(":", 1)[1].strip())
                elif line.startswith("Ref ligand:"):
                    extract_info["Ref_ligand"] = line.split(":", 1)[1].strip()
                elif line.startswith("Scoring function:"):
                    extract_info['scoring_function'] = line.split(":", 1)[1].strip()
                elif line.startswith("output files:"):
                    extract_info['output_file'] = eval(line.split(":", 1)[1].strip())
                elif line.startswith("output logs:"):
                    extract_info['output_log'] = eval(line.split(":", 1)[1].strip())
        
        
        extract_info["extract_receptor_path"] = os.path.normpath(os.path.join(extract_info["extract_dirpath"], f'{extract_info["extract_receptor"]}.pdbqt'))
        extract_info["notfound_receptor_path"] = os.path.normpath(os.path.join(self.current_directory, f'{extract_info["extract_receptor"]}.pdbqt'))
        
        extract_info["extract_refligand_path"] = os.path.normpath(os.path.join(extract_info["extract_dirpath"], f'{extract_info["Ref_ligand"]}.pdbqt'))
        extract_info["notfound_refligand_path"] = os.path.normpath(os.path.join(self.current_directory, f'{extract_info["Ref_ligand"]}.pdbqt'))
        
        #載入Receptor至pymol和表格---------------------------------------------------------------------------------------------------------------
        if os.path.exists(extract_info["extract_receptor_path"]):
            self.load_file_to_pymol(extract_info["extract_receptor_path"])
        elif os.path.exists(extract_info["notfound_receptor_path"]):
            self.load_file_to_pymol(extract_info["notfound_receptor_path"])
        else:
            print("Receptor file was lost, please check your file path.")
        
        if extract_info["Ref_ligand"] != "":
            Ref_ligand_is_uesd = True
            if os.path.exists(extract_info["extract_refligand_path"]):
                self.load_file_to_pymol(extract_info["extract_refligand_path"])
            elif os.path.exists(extract_info["notfound_refligand_path"]):
                self.load_file_to_pymol(extract_info["notfound_refligand_path"])
            else:
                print("Ref ligand file was lost, please check your file path.")
        else:
            Ref_ligand_is_uesd = False
        
        
        self.show_in_receptor_table(extract_info["extract_receptor"], extract_info["Ref_ligand"], Ref_ligand_is_uesd)
        
 
        #載入Ligands至pymol和表格---------------------------------------------------------------------------------------------------------------
        
        for dirname, output_file, log_file in zip(extract_info["extract_ligands"], extract_info['output_file'], extract_info['output_log']):
            each_ligand_pdbqt_path = os.path.normpath(os.path.join(extract_info["extract_dirpath"], dirname, output_file))
            each_log_path = os.path.normpath(os.path.join(extract_info["extract_dirpath"], dirname, log_file))
            if os.path.exists(each_ligand_pdbqt_path):
                result_data = self.extract_affinity_split_pdbqt(dirname, each_ligand_pdbqt_path, each_log_path)
                result_info[dirname] = result_data
                self.load_initial_pdbqt_in_pymol(result_info[dirname])
            elif not os.path.exists(each_ligand_pdbqt_path):
                notfound_ligand_path = os.path.normpath(os.path.join(self.current_directory, dirname, output_file))
                notfound_log_path = os.path.normpath(os.path.join(dirname, self.current_directory, log_file))
                result_data = self.extract_affinity_split_pdbqt(dirname, notfound_ligand_path, notfound_log_path)
                result_info[dirname] = result_data
                self.load_initial_pdbqt_in_pymol(result_info[dirname])
            else:
                print("No such files found, please check your files path are correct.")
        
        
        
        
        
        self.show_in_analysis_ligands_table(extract_info, result_info)
        
        #聚焦回receptor-------------------------------------------------------------------------------------------------------------------------
        if " " in extract_info["extract_receptor"]:
            receptor_name_in_pymol = extract_info["extract_receptor"].replace(" ", "_")
        else:
            receptor_name_in_pymol = extract_info["extract_receptor"]
        self.pymol_process.cmd.zoom(receptor_name_in_pymol)
    
    
    def load_txt_data(self):
        print("load log.txt")
       
        
        
    def data_addition(self):
        print("data added")
    
    
    def show_in_receptor_table(self, receptor_name, reflig_name, reflig_singal):
        receptor_name = receptor_name
        reflig_name = reflig_name
        Ref_ligand_is_uesd = reflig_singal
        
        self.ui.tableWidget_analysis_receptor.setRowCount(len(self.data_dict))  # 確保至少有一行
        

        # 創建一個 QWidget 包含 QCheckBox
        receptor_visible_widget = QWidget()
        receptor_visible_checkbox = QCheckBox()
        receptor_visible_checkbox.setChecked(True)  # 預設選中
        receptor_visible_checkbox.setObjectName("Receptor_Checkbox")  # 設定唯一名稱
        ref_ligand_visible_checkbox = QCheckBox()
        ref_ligand_visible_checkbox.setChecked(True)  # 預設選中
        ref_ligand_visible_checkbox.setObjectName("Ref_Ligand_Checkbox")  # 設定唯一名稱
 
        # 創建一個 QLabel 作為分隔符，顯示斜線 /
        separator = QLabel("/")
        
        # 將 QCheckBox 和 QLabel 添加到布局中
        receptor_visible_layout = QHBoxLayout()
        receptor_visible_layout.addWidget(receptor_visible_checkbox)
        receptor_visible_layout.addWidget(separator)
        receptor_visible_layout.addWidget(ref_ligand_visible_checkbox)
        
        # 調整布局，使控件居中對齊
        receptor_visible_layout.setAlignment(Qt.AlignCenter)  # 居中對齊
        receptor_visible_layout.setContentsMargins(0, 0, 0, 0)  # 設置無邊距
        
        # 將布局應用到 QWidget
        receptor_visible_widget.setLayout(receptor_visible_layout)
        self.ui.tableWidget_analysis_receptor.setCellWidget(0, 2, receptor_visible_widget)
        self.ui.tableWidget_analysis_receptor.setItem(0, 0, QTableWidgetItem(receptor_name))
            
        if Ref_ligand_is_uesd == True:
            ref_ligand_visible_checkbox.setEnabled(True)
            self.ui.tableWidget_analysis_receptor.setItem(0, 1, QTableWidgetItem(reflig_name))
        elif Ref_ligand_is_uesd == False:
            ref_ligand_visible_checkbox.setEnabled(False)
            self.ui.tableWidget_analysis_receptor.setItem(0, 1, QTableWidgetItem("None"))
            
        # 連接 QCheckBox 的信號，當狀態改變時觸發
        receptor_visible_checkbox.stateChanged.connect(lambda: self.visible_signal(receptor_visible_checkbox, receptor_name))
        ref_ligand_visible_checkbox.stateChanged.connect(lambda: self.visible_signal(ref_ligand_visible_checkbox, reflig_name))
    
        
    def show_in_analysis_ligands_table(self, extract_info, result_data_dict):
        extract_info = extract_info
        ligands_list = extract_info["extract_ligands"]
        extract_directory = extract_info["extract_dirpath"]
        result_data_dict = result_data_dict
        
        ligands_amount = len(ligands_list)
        
        self.ui.tableWidget_analysis_ligands.setRowCount(ligands_amount)
        

        # 針對特定列設置不同的調整模式
        header_analysis_ligands = self.ui.tableWidget_analysis_ligands.horizontalHeader()
        header_analysis_ligands.setSectionResizeMode(0, QHeaderView.ResizeToContents) 
        
        
        self.ui.tableWidget_analysis_ligands.itemClicked.connect(lambda item: self.zoom_on_click_result(item, result_data_dict))
        
        
        
        for i, name in enumerate(ligands_list):
            self.ui.tableWidget_analysis_ligands.setItem(i, 0, QTableWidgetItem(name))
            current_ligand_detail = result_data_dict[name]  #當前ligand結果的資料
            
            # 創建自定義的 AffinitySelector 控件，並填入對應的數值
            if name in result_data_dict:
                affinity_widget = AffinitySelector(name, current_ligand_detail)
                
                # 連接信號到主函數 update_pymol_model
                affinity_widget.affinity_changed.connect(lambda ligand_name, mode, result_data_dict: self.update_pymol_model(ligand_name, mode, result_data_dict))

                
                # 設定在表格中
                self.ui.tableWidget_analysis_ligands.setCellWidget(i, 1, affinity_widget)
            
            
            # 為每一行創建一個新的 QWidget 和 QCheckBox
            ligands_visible_widget = QWidget()
            ligands_visible_checkbox = QCheckBox()
            ligands_visible_checkbox.setChecked(True)  # 預設選中
    
            # 將 QCheckBox 添加到布局中
            ligands_visible_layout = QHBoxLayout()
            ligands_visible_layout.addWidget(ligands_visible_checkbox)
    
            # 調整布局，使控件居中對齊
            ligands_visible_layout.setAlignment(Qt.AlignCenter)  # 居中對齊
            ligands_visible_layout.setContentsMargins(0, 0, 0, 0)  # 設置無邊距
    
            # 將布局應用到 QWidget
            ligands_visible_widget.setLayout(ligands_visible_layout)
            
            # 將 QWidget 添加到當前行的單元格中
            self.ui.tableWidget_analysis_ligands.setCellWidget(i, 3, ligands_visible_widget)
            
            # 使用默認參數方式，將變量值傳遞給 lambda，防止變量捕獲問題
            ligands_visible_checkbox.stateChanged.connect(
                lambda state, checkbox=ligands_visible_checkbox, name=name, in_pymol_name=current_ligand_detail["in_pymol_name"] : self.visible_signal(checkbox, in_pymol_name)
                )
            
            # log button創建
            log_button_widget = QWidget()
            open_log_button = QPushButton("\U0001F4C3")
            
            # 創建布局並將 QPushButton 添加到布局中
            log_button_layout = QHBoxLayout()
            log_button_layout.addWidget(open_log_button)
            log_button_layout.setAlignment(Qt.AlignCenter)
            log_button_layout.setContentsMargins(0, 0, 0, 0)
            
            # 創建 QWidget 包含布局，並將其添加到表格第3列
            log_button_widget.setLayout(log_button_layout)
            self.ui.tableWidget_analysis_ligands.setCellWidget(i, 2, log_button_widget)
            
            # 連接 QPushButton 點擊事件
            open_log_button.clicked.connect(lambda _, dir_name=name, log_name = current_ligand_detail["log_file"], dir_path = extract_directory: self.log_button_clicked(dir_name, log_name, dir_path))
    
    
    
    def extract_affinity_split_pdbqt(self, dirname, each_ligand_pdbqt_path, each_log_path):
        ligand_name = dirname
        pdbqt_path = each_ligand_pdbqt_path
        log_path = each_log_path
        
        pdbqt_file_basename = os.path.basename(each_ligand_pdbqt_path).split(".")[0]
        log_file = os.path.basename(log_path)
        
        result_data = {"mode":[], "affinity":[], "pdbqt_split":{}, "pdb_cov_split":{}, "in_pymol_name":pdbqt_file_basename, "log_file":log_file}

        # 讀取Log檔案內容
        with open(log_path, "r") as file:
            content = file.read()
            
        # 使用正則表達式匹配 mode 和 affinity
        log_pattern = re.compile(r"^\s*(\d+)\s+([+-]?\d+\.\d+)\s+.*$", re.MULTILINE)   #^行首, \s空白符號, \d任意一個數字, +一個或多個, *零個或多個, $行尾
        
        for match in log_pattern.finditer(content):
            mode = int(match.group(1))       # mode 轉為整數
            affinity = float(match.group(2))  # affinity 轉為浮點數
            result_data["mode"].append(mode)
            result_data["affinity"].append(affinity)
        
        
        # 讀取 pdbqt 文件內容
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as file:
            pdbqt_content = file.read()
        
        # 使用正則表達式匹配每個 `MODEL` 到 `ENDMDL` 的內容
        pdbqt_pattern = re.compile(r"(MODEL\s+(\d+).*?ENDMDL)", re.DOTALL)
        pdbqt_model_dict = {}
        pdb_model_dict = {}
        
        # 遍歷所有匹配到的 `MODEL` 結構
        for match in pdbqt_pattern.finditer(pdbqt_content):
            pdbqt_model_content = match.group(1)  # 整個模型的內容
            model_index = match.group(2)  
            # 將 `model_index` 作為字典的鍵，`model_content` 作為對應值
            pdbqt_model_dict[model_index] = pdbqt_model_content
            
            # 使用 Open Babel 将 pdbqt 转换为 pdb 格式字符串
            ligand_mode_name = ligand_name + str(model_index)
            pdb_model_content = self.pdbqt_to_pdb(pdbqt_model_content, ligand_mode_name)
            if pdb_model_content:  # 如果转换成功
                pdb_model_dict[model_index] = pdb_model_content
            
        # 將所有的 `model_dict` 加入到 `pdbqt_split` 中
        result_data["pdbqt_split"] = pdbqt_model_dict
        result_data["pdb_cov_split"] = pdb_model_dict

        return result_data
    
     
    def visible_signal(self, checkbox, name):
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
    
    def zoom_on_click_result(self, item, ligand_detail):
        ligand_detail = ligand_detail
        
        # 根據點擊的內容執行 PyMOL 的 zoom 指令
        object_name = item.text()
        in_pymol_name = ligand_detail[object_name]["in_pymol_name"]
        
        # 檢查是否需要將空格替換為底線
        if " " in in_pymol_name:
            in_pymol_name = in_pymol_name.replace(" ", "_")
    
        # 執行 PyMOL 的 zoom 命令
        self.pymol_process.cmd.zoom(in_pymol_name)
    
    
    
    def log_button_clicked(self, dir_name, log_name, dir_path):
        dir_name = dir_name
        log_name = log_name
        dir_path = dir_path
        
        log_path = os.path.normpath(os.path.join(dir_path, dir_name, log_name))
        
        
        # 判斷檔案是否存在
        if os.path.exists(log_path):
            log_path_to_open = log_path
        else:
            log_path_to_open = os.path.normpath(os.path.join(self.current_directory, log_name))
        
        # 如果沒有找到檔案，顯示錯誤訊息
        if not os.path.exists(log_path_to_open):
            print(f"Log file not found: {log_path_to_open}")
            
        # 如果 log_viewer_window 不存在，創建一個新視窗
        if not hasattr(self, 'log_viewer_window') or not self.log_viewer_window:
            self.log_viewer_window = log_viewer.Log_viewer(dir_name, log_path_to_open)
            self.log_viewer_window.raise_()
            self.log_viewer_window.activateWindow()
            self.log_viewer_window.show()
            
            # 視窗關閉後重置 log_viewer_window 為 None
            self.log_viewer_window.finished.connect(self.reset_log_viewer_window)
        else:
            # 如果視窗已經存在，則檢查是否已有相同 dir_name 的 tab
            for index in range(self.log_viewer_window.log_viewer_ui.tabWidget.count()):
                if self.log_viewer_window.log_viewer_ui.tabWidget.tabText(index) == dir_name:
                    # 如果找到相同名稱的 tab，切換到該 tab
                    self.log_viewer_window.log_viewer_ui.tabWidget.setCurrentIndex(index)
                    # 確保視窗顯示在前台
                    self.log_viewer_window.raise_()
                    self.log_viewer_window.activateWindow()
                    return

            # 如果沒有相同的 tab，則添加一個新的 tab
            self.log_viewer_window.add_log(dir_name, log_path_to_open)   
            self.log_viewer_window.raise_()
            self.log_viewer_window.activateWindow()
        
    def reset_log_viewer_window(self):
        """重置 log_viewer_window 狀態"""
        self.log_viewer_window = None
    
    
    
    
    
    
    def affinity_filter_search(self):
        min_value = self.ui.doubleSpinBox_low_affinity.value()
        max_value = self.ui.doubleSpinBox_high_affinity.value()
        
        if min_value > max_value:
            error_message = "Warning message: low value is higher tha high value"
            error_window = QtWidgets.QMessageBox()
            error_window.setIcon(QtWidgets.QMessageBox.Critical)
            error_window.setWindowTitle("Input Error")
            error_window.setInformativeText(error_message)
            error_window.setStandardButtons(QtWidgets.QMessageBox.Ok)
            error_window.exec_()  
            return
        
        # 重設所有行為顯示狀態，確保每次篩選操作都是基於完整行列表進行
        self.affinity_filter_reset()
        
        # 遍歷表格中每一行
        for row in range(self.ui.tableWidget_analysis_ligands.rowCount()):
            # 取得第 1 列（即 "Affinity" 列）中的 AffinitySelector 控件
            affinity_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 1)
            
            # 確認該單元格是否為 AffinitySelector 類型
            if isinstance(affinity_widget, AffinitySelector):
                # 獲取目前選擇的 affinity 值，例如 "1: -8.9164" 格式
                current_affinity_text = affinity_widget.get_current_affinity()
                
                # 解析 affinity 值，提取數值部分
                try:
                    # 將格式 "1: -8.9164" 分解為 mode 和 affinity
                    mode, affinity = current_affinity_text.split(": ")
                    affinity_value = float(affinity)
                    
                    # 取得該行的 Checkbox
                    ligand_checkbox_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 3)
                    if ligand_checkbox_widget:
                        ligand_checkbox = ligand_checkbox_widget.findChild(QCheckBox)
                    
                        # 根據篩選條件設定行可見性 & 勾選狀態
                        if min_value <= affinity_value <= max_value:
                            self.ui.tableWidget_analysis_ligands.setRowHidden(row, False)  # 顯示符合條件的行
                            if ligand_checkbox:
                                ligand_checkbox.setChecked(True)  # ✅ 自動勾選篩選成功的行
                        else:
                            self.ui.tableWidget_analysis_ligands.setRowHidden(row, True)  # 隱藏不符合條件的行
                            if ligand_checkbox:
                                ligand_checkbox.setChecked(False)  # ❌ 自動取消選取不符合條件的行
                except ValueError:
                    print(f"解析 Affinity 值失敗: {current_affinity_text}")
    
    
    def affinity_filter_reset(self):
        row_count = self.ui.tableWidget_analysis_ligands.rowCount()
        for row in range(row_count):
            self.ui.tableWidget_analysis_ligands.showRow(row)


    def load_file_to_pymol(self, filepath):
        if self.pymol_process:
            try:
                self.pymol_process.cmd.load(filepath)
                
            except Exception as e:
                print("Error sending command to PyMOL:", e)
    
    
    def load_initial_pdbqt_in_pymol(self, current_result_data):
        current_result_data = current_result_data
        pdbqt_model_dict = current_result_data["pdbqt_split"]
        pdb_model_dict = current_result_data["pdb_cov_split"]
        in_pymol_name = current_result_data["in_pymol_name"]
        
        
        if self.pymol_process:
            try:
                self.pymol_process.cmd.read_pdbstr(pdb_model_dict["1"], in_pymol_name)
                self.pymol_process.cmd.show("sticks", in_pymol_name)
                self.pymol_process.cmd.hide("spheres", in_pymol_name)   
            except Exception as e:
                print("Error sending command to PyMOL:", e)
                self.pymol_process.cmd.read_pdbstr(pdbqt_model_dict["1"], in_pymol_name)
                self.pymol_process.cmd.show("sticks", in_pymol_name)
                self.pymol_process.cmd.hide("spheres", in_pymol_name)

        
    def update_pymol_model(self, ligand_name, mode, result_data_dict):
        """
        根據 ligand 名稱和選擇的 mode 來更新 PyMOL 中的顯示。
        """
        
        current_ligand_detail = result_data_dict
        pdbqt_model_dict = current_ligand_detail["pdbqt_split"]
        pdb_model_dict = current_ligand_detail["pdb_cov_split"]
        in_pymol_name = current_ligand_detail["in_pymol_name"]
        
        # 檢查 mode 是否在模型字典中
        if str(mode) in pdb_model_dict:

            pdbqt_model_content = pdbqt_model_dict[str(mode)]
            pdb_model_content = pdb_model_dict[str(mode)]
            
            try:
                # 重新載入指定的模型
                self.pymol_process.cmd.read_pdbstr(pdb_model_content, in_pymol_name, state=1, finish=1)
                self.pymol_process.cmd.show("sticks", in_pymol_name)
                self.pymol_process.cmd.hide("spheres", in_pymol_name)  
            except Exception as e:
                print("Error sending command to PyMOL:", e)
                self.pymol_process.cmd.read_pdbstr(pdbqt_model_content, in_pymol_name, state=1, finish=1)
                self.pymol_process.cmd.show("sticks", in_pymol_name)
                self.pymol_process.cmd.hide("spheres", in_pymol_name)   
        
     
    
    def pdbqt_to_pdb(self, pdbqt_model_content, compound_name):
        pdbqt_model_content = pdbqt_model_content
    
        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInFormat("pdbqt")
        ob_conversion.SetOutFormat("pdb")
        
        molecular = openbabel.OBMol()
        ob_conversion.ReadString(molecular, pdbqt_model_content)
        pdb_model_content = ob_conversion.WriteString(molecular)
        # **確保 PDB 內容非空**
        if pdb_model_content:
            # 🔹 使用正則表達式替換 `COMPND` 行
            pdb_model_content = re.sub(r"^(COMPND\s+)(UNNAMED)", lambda match: f"{match.group(1)}{compound_name}", pdb_model_content, flags=re.MULTILINE)   # `match.group(1)` 對應 `COMPND

    
        return pdb_model_content
       
    def save_image_action(self):
        """存取 PyMOL 畫面，使用 QDialog 讓 Checkbox 正確對齊"""

        # 創建對話框
        dialog = QDialog()
        dialog.setWindowTitle("Save Image")
    
        # 訊息標籤
        label = QLabel("Do you want to save the current PyMOL view?", dialog)
        label.setStyleSheet("font-size: 12pt;")  # 設定字體大小為 14pt
    
        # 添加 Checkbox
        ray_checkbox = QCheckBox("Enable ray tracing (higher quality but slower)", dialog)
    
        # 水平佈局讓 Checkbox 置中
        checkbox_layout = QHBoxLayout()
        checkbox_layout.addStretch()  # 左側空白
        checkbox_layout.addWidget(ray_checkbox)  # 添加 Checkbox
        checkbox_layout.addStretch()  # 右側空白
    
        # 按鈕
        button_yes = QPushButton("Yes", dialog)
        button_no = QPushButton("No", dialog)
    
        # 水平佈局讓按鈕置中
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(button_yes)
        button_layout.addWidget(button_no)
        button_layout.addStretch()
    
        # 垂直佈局 (讓元素垂直排列)
        layout = QVBoxLayout()
        layout.addWidget(label, alignment=Qt.AlignCenter)
        layout.addLayout(checkbox_layout)  # 放入 Checkbox 佈局
        layout.addLayout(button_layout)  # 放入按鈕佈局
        
        dialog.setLayout(layout)
    
        # 連接按鈕事件
        button_yes.clicked.connect(dialog.accept)
        button_no.clicked.connect(dialog.reject)
    
        # 顯示對話框
        response = dialog.exec_()
    
        if response == QDialog.Rejected:
            return  # 使用者選擇 No，結束函數
    
        # 開啟 QFileDialog 讓使用者選擇儲存位置
        file_path, _ = QFileDialog.getSaveFileName(
            None, "Save Image", "", "PNG Files (*.png);;JPEG Files (*.jpg);;All Files (*)"
        )
        
        if not file_path:
            return  # 使用者未選擇儲存位置
        
        # 如果選擇使用 ray tracing
        if ray_checkbox.isChecked():
            self.pymol_process.cmd.ray()  # 執行 ray tracing，提高畫質
            self.pymol_process.cmd.png(file_path, dpi=300)
        else:
            self.pymol_process.cmd.draw()  # 嘗試重新繪製
            self.pymol_process.cmd.refresh()  # 強制刷新畫面
            self.pymol_process.cmd.png(file_path, dpi=300)

        if os.path.exists(file_path):
            QMessageBox.information(None, "Save Complete", f"Image saved successfully to:\n{file_path}")
        else:
            QMessageBox.critical(None, "Error", "Failed to save the image. Please try again.")



    def save_complex_action(self):
        # 選擇存檔位置
        file_path, _ = QFileDialog.getSaveFileName(
            None, "Save Complex As", "", "PDB Files (*.pdb);;All Files (*)"
        )
        
        if not file_path:
            return  # 使用者取消存檔
    
    
        # 🔹 取得當前 UI 表格的選擇狀態
        selected_receptor = None
        selected_ref_ligand = None
        selected_ligands = []
        upload_path = list(self.data_dict.keys())[0]  # 取得第一個 key(檔案原始上傳路徑)
        
        # **檢查 Receptor 和 Ref Ligand 的 Checkbox**
        receptor_checkbox_widget = self.ui.tableWidget_analysis_receptor.cellWidget(0, 2)
        if receptor_checkbox_widget:
            checkboxes = receptor_checkbox_widget.findChildren(QCheckBox)
        
            # 初始化變數
            receptor_checkbox = None
            ref_ligand_checkbox = None
        
            # 遍歷所有找到的 QCheckBox
            for checkbox in checkboxes:
                if checkbox.objectName() == "Receptor_Checkbox":
                    receptor_checkbox = checkbox
                elif checkbox.objectName() == "Ref_Ligand_Checkbox":
                    ref_ligand_checkbox = checkbox
                    
            # **檢查是否有選取 Receptor**
            if receptor_checkbox and receptor_checkbox.isChecked():
                selected_receptor = self.ui.tableWidget_analysis_receptor.item(0, 0).text()
            # **檢查是否有選取 Ref Ligand**
            if ref_ligand_checkbox and ref_ligand_checkbox.isChecked():
                selected_ref_ligand = self.ui.tableWidget_analysis_receptor.item(0, 1).text()
    
    
    
        # **檢查 Ligand 表格中的狀態**
        for row in range(self.ui.tableWidget_analysis_ligands.rowCount()):
            ligand_checkbox_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 3)
            if ligand_checkbox_widget:
                ligand_checkbox = ligand_checkbox_widget.findChild(QCheckBox)
                if ligand_checkbox and ligand_checkbox.isChecked():
                    ligand_name = self.ui.tableWidget_analysis_ligands.item(row, 0).text()
    
                    # 取得選擇的 Mode
                    affinity_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 1)
                    if isinstance(affinity_widget, AffinitySelector):   #檢查該 Cell Widget 是否為 AffinitySelector 類型
                        selected_mode = affinity_widget.combo_box.currentIndex() + 1  # Mode 是 1-based index
                        # 從 `self.data_dict` 取得該 Ligand 在當前 Mode 下的 PDB 結構
                        ligand_pdb = self.data_dict[upload_path]["ligand_data"][ligand_name]["pdb_cov_split"].get(str(selected_mode), "")
                        ligand_mode_name = f"{ligand_name}_{selected_mode}"
                        ligand_pdb_with_header = f"HEADER    {ligand_mode_name}\n" + ligand_pdb
                        if ligand_pdb:
                            selected_ligands.append(ligand_pdb_with_header)
    
    
    
        # **確保至少有選擇一個 Receptor 或 Ligand**
        if not selected_receptor and not selected_ligands:
            QMessageBox.warning(None, "No Selection", "No receptor or ligand selected for saving.")
            return
    
        # **讀取選擇的 Receptor PDB**
        receptor_pdb_content = ""
        receptor_header_line = ""
        if selected_receptor:
            for receptor_pdbqt_path in [
                self.data_dict[upload_path]["extract_data"]["extract_receptor_path"],
                self.data_dict[upload_path]["extract_data"]["notfound_receptor_path"]
            ]:
                if os.path.exists(receptor_pdbqt_path):
                    with open(receptor_pdbqt_path, "r", encoding="utf-8") as file:
                        receptor_pdbqt_content = file.read()
                    receptor_pdb_content = self.pdbqt_to_pdb(receptor_pdbqt_content, selected_receptor)
                    receptor_header_line = f"HEADER    {selected_receptor}\n"
                    break  # 找到一個可用的就停止
        
        
        # **讀取選擇的 Ref Ligands PDB**
        ref_ligand_pdb_content = ""
        ref_ligand_header_line = ""
        
        if selected_ref_ligand:
            for ref_ligand_pdbqt_path in [
                self.data_dict[upload_path]["extract_data"]["extract_refligand_path"],
                self.data_dict[upload_path]["extract_data"]["notfound_refligand_path"]
            ]:
                if os.path.exists(ref_ligand_pdbqt_path):
                    with open(ref_ligand_pdbqt_path, "r", encoding="utf-8") as file:
                        ref_ligand_pdbqt_content = file.read()
                    ref_ligand_pdb_content = self.pdbqt_to_pdb(ref_ligand_pdbqt_content, selected_ref_ligand)
                    ref_ligand_header_line = f"HEADER    {selected_ref_ligand}\n"
                    break  # 找到可用的檔案後立即停止
                    
        
       
        
        complex_pdb = receptor_header_line + receptor_pdb_content + "\n" + ref_ligand_header_line + ref_ligand_pdb_content + "\n".join(selected_ligands)
        
        with open(file_path, "w") as pdb_file:
            pdb_file.write(complex_pdb)
    
        QMessageBox.information(None, "Save Complete", f"Complex saved successfully to:\n{file_path}")
             
        
    def save_ligand_action(self):
        """根據選擇的 Ligand 和 Mode 存為 PDB"""

        # 🔹 取得當前 UI 表格的選擇狀態
        selected_ligands = []
        
        upload_path = list(self.data_dict.keys())[0]  # 取得第一個 key(檔案原始上傳路徑)
        
        for row in range(self.ui.tableWidget_analysis_ligands.rowCount()):
            ligand_checkbox_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 3)
            if ligand_checkbox_widget:
                ligand_checkbox = ligand_checkbox_widget.findChild(QCheckBox)
                if ligand_checkbox and ligand_checkbox.isChecked():
                    ligand_name = self.ui.tableWidget_analysis_ligands.item(row, 0).text()
    
                    # 取得選擇的 Mode
                    affinity_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 1)
                    if isinstance(affinity_widget, AffinitySelector):
                        selected_mode = affinity_widget.combo_box.currentIndex() + 1  # Mode 是 1-based index
    
                        # 從 `self.data_dict` 取得該 Ligand 在當前 Mode 下的 PDB 結構
                        ligand_pdb = self.data_dict[upload_path]["ligand_data"][ligand_name]["pdb_cov_split"].get(str(selected_mode), "")
                        ligand_mode_name = f"{ligand_name}_{selected_mode}"
                        ligand_pdb_with_header = f"HEADER    {ligand_mode_name}\n" + ligand_pdb
                        if ligand_pdb:
                            selected_ligands.append(ligand_pdb_with_header)
        
        
    
        # **確保至少有選擇一個 Ligand**
        if not selected_ligands:
            QMessageBox.warning(None, "No Selection", "No ligand selected for saving.")
            return
    
        # **選擇存檔位置**
        file_path, _ = QFileDialog.getSaveFileName(
            None, "Save Ligand As", "", "PDB Files (*.pdb);;All Files (*)"
        )
        
        if not file_path:
            return  # 使用者取消存檔
    
        # **合併選擇的 Ligand PDB 並儲存**
        with open(file_path, "w") as pdb_file:
            all_ligand_content = "\n".join(selected_ligands)
            pdb_file.write(all_ligand_content)
    
        QMessageBox.information(None, "Save Complete", f"Ligand saved successfully to:\n{file_path}")
        
        
        
 #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       
    def toggle_interaction_analysis(self):
        selected_receptor, selected_ligand = self.get_selected_molecules()
        if selected_receptor == None or selected_ligand == None:
            return
        
        # **儲存 PyMOL 當前狀態**
        self.docking_log_session = self.pymol_process.cmd.get_session()

        # **清除並載入 PyMOL**
        self.pymol_process.cmd.reinitialize()
        
        # **使用 read_pdbstr() 讀取 PDB 字串，而不是 load()**
        self.pymol_process.cmd.read_pdbstr(selected_receptor, "Receptor")
        self.pymol_process.cmd.read_pdbstr(selected_ligand, "Ligand")

        # **計算作用力**
        interactions = self.detect_interactions(selected_receptor, selected_ligand)
        if not interactions or all(len(v) == 0 for v in interactions.values()):
            QMessageBox.information(None, "No Interactions", "No interactions detected between the selected receptor and ligand.")
            # **回復 PyMOL 原始狀態**
            self.pymol_process.cmd.reinitialize()
            self.pymol_process.cmd.set_session(self.docking_log_session)
            return
        self.visualize_interaction_in_pymol(interactions)
        
        # **🔹 清空作用力表格**
        self.ui.tableWidget_interaction_analysis.setRowCount(0)  # 先清除舊的行數
        
        # **🔹 切換到作用力分析頁面**
        self.ui.stackedWidget_analysis.setCurrentWidget(self.ui.page_interaction_analysis)

        # **🔹 切換表格成作用力數據**
        self.show_interaction_table(interactions, selected_receptor, selected_ligand)
        self.ui.label_interaction_receptor
        
        
           
    def back_to_analysis_result(self):    
        self.ui.stackedWidget_analysis.setCurrentWidget(self.ui.page_analysis_basic)
        
        # **回復 PyMOL 原始狀態**
        self.pymol_process.cmd.reinitialize()
        self.pymol_process.cmd.set_session(self.docking_log_session)

        
    def get_selected_molecules(self):
        selected_receptor = None
        selected_ref_ligand = None
        selected_ligands = []
        upload_path = list(self.data_dict.keys())[0]  # 取得第一個 key(檔案原始上傳路徑)
        selected_receptor_name = self.ui.tableWidget_analysis_receptor.item(0, 0).text()
        selected_ref_ligand_name = self.ui.tableWidget_analysis_receptor.item(0, 1).text()
        
        # **檢查 Receptor 和 Ref Ligand 的 Checkbox**
        receptor_checkbox_widget = self.ui.tableWidget_analysis_receptor.cellWidget(0, 2)
        if receptor_checkbox_widget:
            checkboxes = receptor_checkbox_widget.findChildren(QCheckBox)
        
            # 初始化變數
            receptor_checkbox = None
            ref_ligand_checkbox = None
        
            # 遍歷所有找到的 QCheckBox
            for checkbox in checkboxes:
                if checkbox.objectName() == "Receptor_Checkbox":
                    receptor_checkbox = checkbox
                elif checkbox.objectName() == "Ref_Ligand_Checkbox":
                    ref_ligand_checkbox = checkbox
                    
            # **檢查是否有選取 Receptor**
            if receptor_checkbox and receptor_checkbox.isChecked(): 
                # **讀取選擇的 Receptor PDB**
                receptor_pdb_content = ""
                receptor_header_line = ""
                for receptor_pdbqt_path in [
                    self.data_dict[upload_path]["extract_data"]["extract_receptor_path"],
                    self.data_dict[upload_path]["extract_data"]["notfound_receptor_path"]
                ]:
                    if os.path.exists(receptor_pdbqt_path):
                        with open(receptor_pdbqt_path, "r", encoding="utf-8") as file:
                            receptor_pdbqt_content = file.read()
                        receptor_pdb_content = self.pdbqt_to_pdb(receptor_pdbqt_content, selected_receptor_name)
                        receptor_header_line = f"HEADER    {selected_receptor_name}\n"
                        selected_receptor = receptor_header_line + receptor_pdb_content
                        break  # 找到一個可用的就停止
  
            # **檢查是否有選取 Ref Ligand**
            if ref_ligand_checkbox and ref_ligand_checkbox.isChecked():
                # **讀取選擇的 Ref Ligands PDB**
                ref_ligand_pdb_content = ""
                ref_ligand_header_line = ""
                for ref_ligand_pdbqt_path in [
                    self.data_dict[upload_path]["extract_data"]["extract_refligand_path"],
                    self.data_dict[upload_path]["extract_data"]["notfound_refligand_path"]
                ]:
                    if os.path.exists(ref_ligand_pdbqt_path):
                        with open(ref_ligand_pdbqt_path, "r", encoding="utf-8") as file:
                            ref_ligand_pdbqt_content = file.read()
                        ref_ligand_pdb_content = self.pdbqt_to_pdb(ref_ligand_pdbqt_content, selected_ref_ligand_name)
                        ref_ligand_header_line = f"HEADER    {selected_ref_ligand_name}\n"
                        selected_ref_ligand = ref_ligand_header_line + ref_ligand_pdb_content
                        break  # 找到可用的檔案後立即停止


        # **檢查 Ligand 表格中的狀態**
        for row in range(self.ui.tableWidget_analysis_ligands.rowCount()):
            ligand_checkbox_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 3)
            if ligand_checkbox_widget:
                ligand_checkbox = ligand_checkbox_widget.findChild(QCheckBox)
                if ligand_checkbox and ligand_checkbox.isChecked():
                    ligand_name = self.ui.tableWidget_analysis_ligands.item(row, 0).text()
    
                    # 取得選擇的 Mode
                    affinity_widget = self.ui.tableWidget_analysis_ligands.cellWidget(row, 1)
                    if isinstance(affinity_widget, AffinitySelector):   #檢查該 Cell Widget 是否為 AffinitySelector 類型
                        selected_mode = affinity_widget.combo_box.currentIndex() + 1  # Mode 是 1-based index
                        # 從 `self.data_dict` 取得該 Ligand 在當前 Mode 下的 PDB 結構
                        ligand_pdb = self.data_dict[upload_path]["ligand_data"][ligand_name]["pdb_cov_split"].get(str(selected_mode), "")
                        ligand_mode_name = f"{ligand_name}_{selected_mode}"
                        ligand_pdb_with_header = f"HEADER    {ligand_mode_name}\n" + ligand_pdb
                        if ligand_pdb:
                            selected_ligands.append(ligand_pdb_with_header)
                            
       
    
        # **確保至少有選擇一個 Receptor 或 Ligand**
        if not selected_receptor:
            QMessageBox.warning(None, "No Selection", "No receptor is picking.")
            return None, None
        elif len(selected_ligands) > 1:
            QMessageBox.warning(None, "Selection invalid", "One ligand only.")
            return None, None
        elif not selected_ligands and not selected_ref_ligand:
            QMessageBox.warning(None, "No Selection", "No ligand is picking.")
            return None, None
        elif selected_ligands and selected_ref_ligand:
            QMessageBox.warning(None, "Selection invalid", "One ligand only.")
            return None, None
        elif not selected_ligands and selected_ref_ligand:
            return selected_receptor, selected_ref_ligand
        elif selected_ligands and not selected_ref_ligand:
            return selected_receptor, str(selected_ligands[0])
        
    
    
    def detect_interactions(self, selected_receptor, selected_ligand):
        """ 計算 Receptor-Ligand 之間的不同類型的相互作用，並輸出 Residue 相關資訊 """
        
        # 解析 PDB，獲取 Residue 資訊
        receptor_residue_map = self.parse_pdb_residues(selected_receptor)
        ligand_residue_map = self.parse_pdb_residues(selected_ligand)
        
        
        
        # **檢查資料是否完整**
        if not receptor_residue_map or not ligand_residue_map:
            print("⚠️ 錯誤：Receptor 或 Ligand 檔案格式不正確")
            return {}  # 直接回傳空字典，避免後續出錯
        
        receptor = Chem.MolFromPDBBlock(selected_receptor, removeHs=False)
        ligand = Chem.MolFromPDBBlock(selected_ligand, removeHs=False)
    
        if receptor is None or ligand is None:
            print("讀取 PDB 失敗")
            return {}
    
        hbond_list = []
        hydrophobic_list = []
        pi_stacking_list = []
        salt_bridge_list = []
    
        # **受體 (Receptor) 作為供氫者 (Donor)**
        for rec_atom in receptor.GetAtoms():
            if rec_atom.GetSymbol() in ["O", "N"]:  
                for neighbor in rec_atom.GetNeighbors():
                    if neighbor.GetSymbol() == "H":  
                        rec_donor_hydrogen = neighbor.GetIdx()
    
                        try:
                            rec_donor_hydrogen_pos = receptor.GetConformer().GetAtomPosition(rec_donor_hydrogen)
                            rec_donor_pos = receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx())
                        except:
                            print(f"錯誤：無法獲取 Receptor 原子 {rec_donor_hydrogen} 或 {rec_atom.GetIdx()} 的座標！")
                            continue  # 跳過錯誤的原子
    
                        for lig_atom in ligand.GetAtoms():
                            if lig_atom.GetSymbol() in ["O", "N"]:
                                try:
                                    lig_acceptor_pos = ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())
                                except:
                                    print(f"錯誤：無法獲取 Ligand 原子 {lig_atom.GetIdx()} 的座標！")
                                    continue  
    
                                h_a_distance = (rec_donor_hydrogen_pos - lig_acceptor_pos).Length()
                                d_a_distance = (rec_donor_pos - lig_acceptor_pos).Length()
    
                                if h_a_distance < 2.5 and d_a_distance < 3.5:
                                    # **透過 XYZ 座標回推 ATOM ID**
                                    xyz_key_rec = (round(rec_donor_pos.x, 3), round(rec_donor_pos.y, 3), round(rec_donor_pos.z, 3))
                                    xyz_key_lig = (round(lig_acceptor_pos.x, 3), round(lig_acceptor_pos.y, 3), round(lig_acceptor_pos.z, 3))
    
                                    rec_header = next(iter(receptor_residue_map))  
                                    rec_atom_id = receptor_residue_map[rec_header]["xyz_to_atom_id"].get(xyz_key_rec, None)
    
                                    lig_header = next(iter(ligand_residue_map))
                                    lig_atom_id = ligand_residue_map[lig_header]["xyz_to_atom_id"].get(xyz_key_lig, None)
    
                                    if rec_atom_id and lig_atom_id:
                                        # 直接透過 ATOM ID 查找完整氨基酸資訊
                                        rec_res_info = receptor_residue_map[rec_header]["atom_id_map"].get(rec_atom_id, {})
                                        lig_res_info = ligand_residue_map[lig_header]["atom_id_map"].get(lig_atom_id, {})
                                        
                                        rec_donor_atom_name = rec_res_info.get('Atom Name', 'UNK atom')
                                        rec_donor_atom_id = rec_res_info.get('ATOM ID', 'UNK id')
                                        lig_accep_atom_name = lig_res_info.get('Atom Name', 'UNK atom')
                                        lig_accep_atom_id = lig_res_info.get('ATOM ID', 'UNK id')

                                        rec_d_to_lig_a_atom = f"R:{rec_donor_atom_name}({rec_donor_atom_id}) → L:{lig_accep_atom_name}({lig_accep_atom_id})"

                                        receptor_residue_name = rec_res_info.get('Residue Name', 'UNK')
                                        receptor_residue_id = rec_res_info.get('Residue ID', 'UNK id')
                                        
                                        receptor_residue_name_id = f"{receptor_residue_name}({receptor_residue_id})"
                                        
                                        distance = f"H--{round(h_a_distance, 2)}--A"

                                        hbond_list.append((rec_d_to_lig_a_atom, receptor_residue_name_id, distance))
    
    
        # **配體 (Ligand) 作為供氫者 (Donor)**
        for lig_atom in ligand.GetAtoms():
            if lig_atom.GetSymbol() in ["O", "N"]:  
                for neighbor in lig_atom.GetNeighbors():
                    if neighbor.GetSymbol() == "H":  
                        lig_donor_hydrogen = neighbor.GetIdx()
    
                        try:
                            lig_donor_hydrogen_pos = ligand.GetConformer().GetAtomPosition(lig_donor_hydrogen)
                            lig_donor_pos = ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())
                        except:
                            print(f"錯誤：無法獲取 Ligand 原子 {lig_donor_hydrogen} 或 {lig_atom.GetIdx()} 的座標！")
                            continue  # 跳過錯誤的原子
    
                        for rec_atom in receptor.GetAtoms():
                            if rec_atom.GetSymbol() in ["O", "N"]:
                                try:
                                    rec_acceptor_pos = receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx())
                                except:
                                    print(f"錯誤：無法獲取 receptor 原子 {rec_atom.GetIdx()} 的座標！")
                                    continue  
    
                                h_a_distance = (lig_donor_hydrogen_pos - rec_acceptor_pos).Length()
                                d_a_distance = (lig_donor_pos - rec_acceptor_pos).Length()
    
                                if h_a_distance < 2.5 and d_a_distance < 3.5:
                                    # **透過 XYZ 座標回推 ATOM ID**
                                    xyz_key_lig = (round(lig_donor_pos.x, 3), round(lig_donor_pos.y, 3), round(lig_donor_pos.z, 3))
                                    xyz_key_rec = (round(rec_acceptor_pos.x, 3), round(rec_acceptor_pos.y, 3), round(rec_acceptor_pos.z, 3))
    
                                    rec_header = next(iter(receptor_residue_map))  
                                    rec_atom_id = receptor_residue_map[rec_header]["xyz_to_atom_id"].get(xyz_key_rec, None)
    
                                    lig_header = next(iter(ligand_residue_map))
                                    lig_atom_id = ligand_residue_map[lig_header]["xyz_to_atom_id"].get(xyz_key_lig, None)
    
                                    if rec_atom_id and lig_atom_id:
                                        # 直接透過 ATOM ID 查找完整氨基酸資訊
                                        rec_res_info = receptor_residue_map[rec_header]["atom_id_map"].get(rec_atom_id, {})
                                        lig_res_info = ligand_residue_map[lig_header]["atom_id_map"].get(lig_atom_id, {})
                                        
                                        lig_donor_atom_name = lig_res_info.get('Atom Name', 'UNK atom')
                                        lig_donor_atom_id = lig_res_info.get('ATOM ID', 'UNK id')
                                        rec_accept_atom_name = rec_res_info.get('Atom Name', 'UNK atom')
                                        rec_accept_atom_id = rec_res_info.get('ATOM ID', 'UNK id')

                                        lig_d_to_rec_a_atom = f"L:{lig_donor_atom_name}({lig_donor_atom_id}) → R:{rec_accept_atom_name}({rec_accept_atom_id})"
                                        
                                        receptor_residue_name = rec_res_info.get('Residue Name', 'UNK')
                                        receptor_residue_id = rec_res_info.get('Residue ID', 'UNK id')
                                        
                                        receptor_residue_name_id = f"{receptor_residue_name}({receptor_residue_id})"
                                        
                                        distance = f"H--{round(h_a_distance, 2)}--A"

                                        hbond_list.append((lig_d_to_rec_a_atom, receptor_residue_name_id, distance))
                                        
                                        
        # **疏水作用**
        for rec_atom in receptor.GetAtoms():
            for lig_atom in ligand.GetAtoms():
                distance = (receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx()) -
                            ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())).Length()

                if distance < 5.0 and rec_atom.GetSymbol() == "C" and lig_atom.GetSymbol() == "C":
                    rec_atom_pos = receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx())
                    lig_atom_pos = ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())
                    
                    xyz_key_rec = (round(rec_atom_pos.x, 3), round(rec_atom_pos.y, 3), round(rec_atom_pos.z, 3))
                    xyz_key_lig = (round(lig_atom_pos.x, 3), round(lig_atom_pos.y, 3), round(lig_atom_pos.z, 3))
                    
                    rec_header = next(iter(receptor_residue_map))  
                    rec_atom_id = receptor_residue_map[rec_header]["xyz_to_atom_id"].get(xyz_key_rec, None)

                    lig_header = next(iter(ligand_residue_map))
                    lig_atom_id = ligand_residue_map[lig_header]["xyz_to_atom_id"].get(xyz_key_lig, None)
                    
                    if rec_atom_id is None or lig_atom_id is None:  
                        continue
                    
                    rec_res_info = receptor_residue_map[rec_header]["atom_id_map"].get(rec_atom_id, {})
                    lig_res_info = ligand_residue_map[lig_header]["atom_id_map"].get(lig_atom_id, {})
                    
                    rec_hydrophobic_atom_name = rec_res_info.get('Atom Name', 'UNK atom')
                    rec_hydrophobic_atom_id = rec_res_info.get('ATOM ID', 'UNK id')
                    lig_hydrophobic_atom_name = lig_res_info.get('Atom Name', 'UNK atom')
                    lig_hydrophobic_atom_id = lig_res_info.get('ATOM ID', 'UNK id')
                    
                    receptor_residue_name = rec_res_info.get('Residue Name', 'UNK')
                    receptor_residue_id = rec_res_info.get('Residue ID', 'UNK id')
                    
                    rec_lig_atom_connect = f"{rec_hydrophobic_atom_name}({rec_hydrophobic_atom_id}) ↔ {lig_hydrophobic_atom_name}({lig_hydrophobic_atom_id})"
                    receptor_residue_name_id = f"{receptor_residue_name}({receptor_residue_id})"

                    hydrophobic_list.append((rec_lig_atom_connect, receptor_residue_name_id, round(distance, 2)))

        # **Pi-Stacking**
        for rec_atom in receptor.GetAtoms():
            for lig_atom in ligand.GetAtoms():
                distance = (receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx()) -
                            ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())).Length()

                if distance < 6.0 and rec_atom.GetIsAromatic() and lig_atom.GetIsAromatic():
                    rec_atom_pos = receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx())
                    lig_atom_pos = ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())
                    
                    xyz_key_rec = (round(rec_atom_pos.x, 3), round(rec_atom_pos.y, 3), round(rec_atom_pos.z, 3))
                    xyz_key_lig = (round(lig_atom_pos.x, 3), round(lig_atom_pos.y, 3), round(lig_atom_pos.z, 3))
                    
                    rec_header = next(iter(receptor_residue_map))  
                    rec_atom_id = receptor_residue_map[rec_header]["xyz_to_atom_id"].get(xyz_key_rec, None)

                    lig_header = next(iter(ligand_residue_map))
                    lig_atom_id = ligand_residue_map[lig_header]["xyz_to_atom_id"].get(xyz_key_lig, None)
                    
                    if rec_atom_id is None or lig_atom_id is None:  
                        continue
                    
                    rec_res_info = receptor_residue_map[rec_header]["atom_id_map"].get(rec_atom_id, {})
                    lig_res_info = ligand_residue_map[lig_header]["atom_id_map"].get(lig_atom_id, {})
                    
                    rec_pistack_atom_name = rec_res_info.get('Atom Name', 'UNK atom')
                    rec_pistack_atom_id = rec_res_info.get('ATOM ID', 'UNK id')
                    lig_pistack_atom_name = lig_res_info.get('Atom Name', 'UNK atom')
                    lig_pistack_atom_id = lig_res_info.get('ATOM ID', 'UNK id')
                    
                    receptor_residue_name = rec_res_info.get('Residue Name', 'UNK')
                    receptor_residue_id = rec_res_info.get('Residue ID', 'UNK id')
                    
                    rec_lig_atom_connect = f"{rec_pistack_atom_name}({rec_pistack_atom_id}) ↔ {lig_pistack_atom_name}({lig_pistack_atom_id})"
                    receptor_residue_name_id = f"{receptor_residue_name}({receptor_residue_id})"

                    pi_stacking_list.append((rec_lig_atom_connect, receptor_residue_name_id, round(distance, 2)))

        # **鹽橋**
        for rec_atom in receptor.GetAtoms():
            for lig_atom in ligand.GetAtoms():
                distance = (receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx()) -
                            ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())).Length()

                if distance < 4.0 and abs(rec_atom.GetFormalCharge()) > 0 and abs(lig_atom.GetFormalCharge()) > 0:
                    rec_atom_pos = receptor.GetConformer().GetAtomPosition(rec_atom.GetIdx())
                    lig_atom_pos = ligand.GetConformer().GetAtomPosition(lig_atom.GetIdx())
                    
                    xyz_key_rec = (round(rec_atom_pos.x, 3), round(rec_atom_pos.y, 3), round(rec_atom_pos.z, 3))
                    xyz_key_lig = (round(lig_atom_pos.x, 3), round(lig_atom_pos.y, 3), round(lig_atom_pos.z, 3))
                    
                    rec_header = next(iter(receptor_residue_map))  
                    rec_atom_id = receptor_residue_map[rec_header]["xyz_to_atom_id"].get(xyz_key_rec, None)

                    lig_header = next(iter(ligand_residue_map))
                    lig_atom_id = ligand_residue_map[lig_header]["xyz_to_atom_id"].get(xyz_key_lig, None)
                    
                    if rec_atom_id is None or lig_atom_id is None:  
                        continue
                    
                    rec_res_info = receptor_residue_map[rec_header]["atom_id_map"].get(rec_atom_id, {})
                    lig_res_info = ligand_residue_map[lig_header]["atom_id_map"].get(lig_atom_id, {})
                    
                    rec_saltbridge_atom_name = rec_res_info.get('Atom Name', 'UNK atom')
                    rec_saltbridge_atom_id = rec_res_info.get('ATOM ID', 'UNK id')
                    lig_saltbridge_atom_name = lig_res_info.get('Atom Name', 'UNK atom')
                    lig_saltbridge_atom_id = lig_res_info.get('ATOM ID', 'UNK id')
                    
                    receptor_residue_name = rec_res_info.get('Residue Name', 'UNK')
                    receptor_residue_id = rec_res_info.get('Residue ID', 'UNK id')
                    
                    rec_lig_atom_connect = f"{rec_saltbridge_atom_name}({rec_saltbridge_atom_id}) ↔ {lig_saltbridge_atom_name}({lig_saltbridge_atom_id})"
                    receptor_residue_name_id = f"{receptor_residue_name}({receptor_residue_id})"

                    salt_bridge_list.append((rec_lig_atom_connect, receptor_residue_name_id, round(distance, 2)))      
                    
                    

        interactions = {
            "H-Bond": hbond_list,
            "Hydrophobic": hydrophobic_list,
            "Pi-Stacking": pi_stacking_list,
            "Salt Bridge": salt_bridge_list
        }                                
        
        
        
        
       
        

        return interactions
    
    def parse_pdb_residues(self, pdb_text):
        pdb_map = {}
        header_regex = re.compile(r"^HEADER\s+(.+)$")
    
        pdb_atom_regex = re.compile(
            r"^(ATOM|HETATM)\s+(\d+)\s+(\S+)\s+(\S+)\s*(\S?)\s*(\d*)\s+"  # ATOM ID, Atom Name, Residue Name, Chain ID (可選), Residue ID
            r"(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+"  # X, Y, Z
            r"(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\S+)\s*$"  # Occupancy, B-Factor, Element
        )
    
        header_key = "UNNAMED"  # 預設 key，防止 PDB 檔沒有 HEADER
        for line in pdb_text.split("\n"):
            header_match = header_regex.match(line)
            info_match = pdb_atom_regex.match(line)
    
            if header_match:
                header_key = header_match.group(1).strip()  # 取 HEADER 內容
                if header_key not in pdb_map:
                    pdb_map[header_key] = {"xyz_to_atom_id": {}, "atom_id_map": {}}
    
            # **解析 ATOM / HETATM 行**
            if info_match:
                atom_id = int(info_match.group(2))  # ATOM ID 總是數字
                atom_name = info_match.group(3).strip()
                residue_name = info_match.group(4).strip()
    
                # **修正 Chain ID (防止 None 或 空白)**
                chain_id = info_match.group(5)
                if chain_id is None or chain_id.strip() == "":
                    chain_id = "UNK"  # 設定未知鏈標識符
    
                # **修正 Residue ID**
                residue_id = info_match.group(6).strip()
                if residue_id.isdigit():
                    residue_id = int(residue_id)  # 轉換為數字
                else:
                    residue_id = "UNK"  # 設為 UNK
    
                # **坐標資訊**
                x, y, z = float(info_match.group(7)), float(info_match.group(8)), float(info_match.group(9))
                element = info_match.group(12).strip()
    
                atom_data = {
                    "ATOM ID": atom_id,
                    "Atom Name": atom_name,
                    "Residue Name": residue_name,
                    "Chain ID": chain_id,
                    "Residue ID": residue_id,
                    "X": x,
                    "Y": y,
                    "Z": z,
                    "Element Symbol": element,
                }
    
                # **XYZ → ATOM ID 映射**
                xyz_key = (x, y, z)
                pdb_map[header_key]["xyz_to_atom_id"][xyz_key] = (chain_id, atom_id)  # 加入鏈標識符，防止重複
    
                # **ATOM ID → 詳細資訊映射**
                atom_key = (chain_id, atom_id)  # 使用 (Chain ID, ATOM ID) 作為 Key，確保不同鏈的相同 ATOM ID 不會覆蓋
                pdb_map[header_key]["atom_id_map"][atom_key] = atom_data
    
        return pdb_map
        
        
        
    def show_interaction_table(self, interactions, selected_receptor, selected_ligand):
        """顯示 Interaction Analysis 的結果表格"""
        for line in selected_receptor.splitlines():
            line = line.strip()  # 去掉前後空白
            if line.startswith("HEADER"):
                parts = line.split(maxsplit=1)  # 只分割一次
                if len(parts) > 1:
                    header_line = parts[1].strip()  # 取得 HEADER 內容
                    self.ui.label_interaction_receptor.setText(header_line)
            break
        
        for line in selected_ligand.splitlines():
            line = line.strip()  # 去掉前後空白
            if line.startswith("HEADER"):
                parts = line.split(maxsplit=1)  # 只分割一次
                if len(parts) > 1:
                    header_line = parts[1].strip()  # 取得 HEADER 內容
                    self.ui.label_interaction_ligand.setText(header_line)
            break
  
        
        for interaction_type, bonds in interactions.items():
            for bond in bonds:
                row = self.ui.tableWidget_interaction_analysis.rowCount()
                self.ui.tableWidget_interaction_analysis.insertRow(row)
                
                atom, residue, distance = bond  
                self.ui.tableWidget_interaction_analysis.setItem(row, 0, QTableWidgetItem(interaction_type))
                self.ui.tableWidget_interaction_analysis.setItem(row, 1, QTableWidgetItem(str(atom)))
                self.ui.tableWidget_interaction_analysis.setItem(row, 2, QTableWidgetItem(str(residue)))  
                self.ui.tableWidget_interaction_analysis.setItem(row, 3, QTableWidgetItem(str(distance)))
    
 
    
                # **🔹 加入 CheckBox 控制作用力顯示**
                checkbox_widget = QWidget()
                checkbox = QCheckBox()
                # **預設 CheckBox 狀態**
                if interaction_type == "H-Bond":
                    checkbox.setChecked(True)   # H-Bond 預設顯示
                else:
                    checkbox.setChecked(False)  # 其他作用力預設關閉
    
                layout = QHBoxLayout()
                layout.addWidget(checkbox)
                layout.setAlignment(Qt.AlignCenter)
                layout.setContentsMargins(0, 0, 0, 0)
                checkbox_widget.setLayout(layout)
    
                self.ui.tableWidget_interaction_analysis.setCellWidget(row, 4, checkbox_widget)
                
                
                if interaction_type == "H-Bond":
                    pattern = re.compile(r"[LR]:\w+\((\d+)\) → [LR]:\w+\((\d+)\)")
                    match = pattern.match(atom)
                    donor_atom_id = match.group(1).strip()
                    acceptor_atom_id = match.group(2).strip()
                    
                    pymol_name = f"{interaction_type}_{donor_atom_id}_{acceptor_atom_id}"
                
                elif interaction_type == "Hydrophobic":
                    pattern = re.compile(r"(\w+)\((\d+)\) ↔ (\w+)\((\d+)\)")
                    match = pattern.match(atom)
                    receptor_hydrophobic_atom_id = match.group(2).strip()
                    ligand_hydrophobic_atom_id = match.group(4).strip()
                    
                    pymol_name = f"{interaction_type}_{receptor_hydrophobic_atom_id}_{ligand_hydrophobic_atom_id}"
                    
                elif interaction_type == "Pi-Stacking":
                    pattern = re.compile(r"(\w+)\((\d+)\) ↔ (\w+)\((\d+)\)")
                    match = pattern.match(atom)
                    receptor_pistack_atom_id = match.group(2).strip()
                    ligand_pistack_atom_id = match.group(4).strip()
                    
                    pymol_name = f"{interaction_type}_{receptor_pistack_atom_id}_{ligand_pistack_atom_id}"
                    
                elif interaction_type == "Salt Bridge":
                    pattern = re.compile(r"(\w+)\((\d+)\) ↔ (\w+)\((\d+)\)")
                    match = pattern.match(atom)
                    receptor_saltbridge_atom_id = match.group(2).strip()
                    ligand_saltbridge_atom_id = match.group(4).strip()
                
                    pymol_name = f"{interaction_type}_{receptor_saltbridge_atom_id}_{ligand_saltbridge_atom_id}"
                
                
                # **🔹 連結 CheckBox 控制 PyMOL 顯示**
                checkbox.stateChanged.connect(
                    lambda state, name=pymol_name: self.toggle_pymol_interaction(state, name)
                )
        
        
        


    def visualize_interaction_in_pymol(self, interactions):
        """ 讓 PyMOL 可視化 H-Bond 作用力 """
        colors = {
            "H-Bond": "yellow",
            "Hydrophobic": "orange",
            "Pi-Stacking": "blue",
            "Salt Bridge": "purple",
        }
        
        pymol_visualize_name = []   #初始不顯示放進這個列表
        for interaction_type, bonds in interactions.items():
            if interaction_type == "H-Bond":  # 只處理氫鍵
                for atom_connect_info, rec_residue_info, distance in bonds:
                    pattern = re.compile(r"[LR]:\w+\((\d+)\) → [LR]:\w+\((\d+)\)")
                    match = pattern.match(atom_connect_info)
                    donor_atom_id = match.group(1).strip()
                    acceptor_atom_id = match.group(2).strip()
                    
                    pymol_name = f"{interaction_type}_{donor_atom_id}_{acceptor_atom_id}"
                    
                    if atom_connect_info.startswith("R:"):
                        self.pymol_process.cmd.distance(
                            pymol_name, 
                            f"Receptor and id {donor_atom_id}", 
                            f"Ligand and id {acceptor_atom_id}"
                        )
                        
                    else:
                        self.pymol_process.cmd.distance(
                            pymol_name, 
                            f"Receptor and id {acceptor_atom_id}", 
                            f"Ligand and id {donor_atom_id}"
                        )
                    
                    
                    self.pymol_process.cmd.show("lines", f"Receptor and resi {rec_residue_info.split('(')[-1].strip(')')}") 
                    self.pymol_process.cmd.set("dash_color", colors[interaction_type], pymol_name)
                    
            elif interaction_type == "Hydrophobic":
                for atom_connect_info, rec_residue_info, distance in bonds:
                    pattern = re.compile(r"(\w+)\((\d+)\) ↔ (\w+)\((\d+)\)")
                    match = pattern.match(atom_connect_info)
                    receptor_hydrophobic_atom_id = match.group(2).strip()
                    ligand_hydrophobic_atom_id = match.group(4).strip()
                    
                    pymol_name = f"{interaction_type}_{receptor_hydrophobic_atom_id}_{ligand_hydrophobic_atom_id}"
                    
                    self.pymol_process.cmd.distance(
                        pymol_name, 
                        f"Receptor and id {receptor_hydrophobic_atom_id}", 
                        f"Ligand and id {ligand_hydrophobic_atom_id}"
                    )
                    self.pymol_process.cmd.set("dash_color", colors[interaction_type], pymol_name)
                    pymol_visualize_name.append(pymol_name)  # 儲存這些物件名稱
                
            elif interaction_type == "Pi-Stacking":
                for atom_connect_info, rec_residue_info, distance in bonds:
                    pattern = re.compile(r"(\w+)\((\d+)\) ↔ (\w+)\((\d+)\)")
                    match = pattern.match(atom_connect_info)
                    receptor_pistack_atom_id = match.group(2).strip()
                    ligand_pistack_atom_id = match.group(4).strip()
                    
                    pymol_name = f"{interaction_type}_{receptor_pistack_atom_id}_{ligand_pistack_atom_id}"
                    
                    self.pymol_process.cmd.distance(
                        pymol_name, 
                        f"Receptor and id {receptor_pistack_atom_id}", 
                        f"Ligand and id {ligand_pistack_atom_id}"
                    )
                    self.pymol_process.cmd.set("dash_color", colors[interaction_type], pymol_name)
                    pymol_visualize_name.append(pymol_name)  # 儲存這些物件名稱
                
            elif interaction_type == "Salt Bridge":
                for atom_connect_info, rec_residue_info, distance in bonds:
                    pattern = re.compile(r"(\w+)\((\d+)\) ↔ (\w+)\((\d+)\)")
                    match = pattern.match(atom_connect_info)
                    receptor_saltbridge_atom_id = match.group(2).strip()
                    ligand_saltbridge_atom_id = match.group(4).strip()
                
                    pymol_name = f"{interaction_type}_{receptor_saltbridge_atom_id}_{ligand_saltbridge_atom_id}"
                    
                    self.pymol_process.cmd.distance(
                        pymol_name, 
                        f"Receptor and id {receptor_saltbridge_atom_id}", 
                        f"Ligand and id {ligand_saltbridge_atom_id}"
                    )
                    self.pymol_process.cmd.set("dash_color", colors[interaction_type], pymol_name)
                    pymol_visualize_name.append(pymol_name)  # 儲存這些物件名稱
    
        # **關閉所有剛載入的作用力物件**
        for pymol_name in pymol_visualize_name:
            self.pymol_process.cmd.disable(pymol_name)   
            
        self.pymol_process.cmd.zoom("Ligand")
    
    
    def toggle_pymol_interaction(self, state, pymol_name):
        """ 控制 PyMOL 作用力顯示 """
        if state == Qt.Checked:
            self.pymol_process.cmd.enable(pymol_name)  # 顯示
        else:
            self.pymol_process.cmd.disable(pymol_name)  # 隱藏
    
    
    
    def save_interaction_action(self):
        """儲存完整的 PyMOL session，包括所有物件"""
        file_path, _ = QFileDialog.getSaveFileName(
            None, "Save PyMOL Session", "", "PyMOL Session (*.pse);;All Files (*)"
        )
    
        if not file_path:
            return  # 使用者取消存檔
    
        # 確保 PyMOL 物件存在
        if not hasattr(self.pymol_process, "cmd"):
            QMessageBox.warning(None, "PyMOL Error", "PyMOL process is not initialized.")
            return
    
        # **儲存 PyMOL Session**
        try:
            self.pymol_process.cmd.save(file_path)
            QMessageBox.information(None, "Save Complete", f"Full session saved:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(None, "Save Error", f"Failed to save session:\n{e}")
            
        
    
    
    
    
    def header_clicked(self, index, table):
        if table == self.ui.tableWidget_analysis_receptor:
            self.ana_receptor_header_vis_state = not self.ana_receptor_header_vis_state
            new_state = self.ana_receptor_header_vis_state
            # 更改表頭的圖標
            self.update_header_icon(table, new_state)

            # 遍歷所有行的 Checkbox，並設置其狀態
            for row in range(table.rowCount()):
                checkbox_widget = table.cellWidget(row, 2)  # 假設 Checkbox 在第 1 列
                if checkbox_widget:
                    receptor_checkbox = checkbox_widget.findChildren(QCheckBox)
                    if len(receptor_checkbox) == 2:  # 确保找到两个复选框
                        receptor_checkbox[0].setChecked(not new_state)  # 切换受体 Checkbox 状态
                        receptor_checkbox[1].setChecked(not new_state)  # 切换参考配体 Checkbox 状态
                        
        elif table == self.ui.tableWidget_analysis_ligands:
            self.ana_ligands_header_vis_state = not self.ana_ligands_header_vis_state
            new_state = self.ana_ligands_header_vis_state
            # 更改表頭的圖標
            self.update_header_icon(table, new_state)

            # 遍歷所有行的 Checkbox，並設置其狀態
            for row in range(table.rowCount()):
                checkbox_widget = table.cellWidget(row, 3)  # 假設 Checkbox 在第 1 列
                if checkbox_widget:
                    checkbox = checkbox_widget.findChild(QCheckBox)
                    if checkbox:
                        checkbox.setChecked(not new_state)  # 切換狀態 (開關顯示)
                        
        elif table == self.ui.tableWidget_interaction_analysis:
            self.interaction_analysis_vis_state = not self.interaction_analysis_vis_state
            new_state = self.interaction_analysis_vis_state
            self.update_header_icon(table, new_state)

            if index == 4:
                # 遍歷所有行的 Checkbox，並設置其狀態
                for row in range(table.rowCount()):
                    checkbox_widget = table.cellWidget(row, 4)  
                    if checkbox_widget:
                        checkbox = checkbox_widget.findChild(QCheckBox)
                        if checkbox:
                            checkbox.setChecked(not new_state)  # 切換狀態 (開關顯示)
                
        
            
                            
    def update_header_icon(self, table, state):
        header_item = QTableWidgetItem()
        if table == self.ui.tableWidget_interaction_analysis:
            if state == False:
                header_item.setText("👁️")  # 當顯示時，設置表頭為「眼睛」圖案
            elif state == True:
                header_item.setText("︶")  # 當隱藏時，設置表頭為「隱藏」圖案
    
            # 假設你要更新的是第二列（第1列，因為索引從0開始）
            table.setHorizontalHeaderItem(4, header_item)
            
        else:    
            if state == False:
                header_item.setText("👁️")  # 當顯示時，設置表頭為「眼睛」圖案
            elif state == True:
                header_item.setText("︶")  # 當隱藏時，設置表頭為「隱藏」圖案
    
            # 假設你要更新的是第二列（第1列，因為索引從0開始）
            table.setHorizontalHeaderItem(3, header_item)
            
            
    

    
    
class AffinitySelector(QWidget):
    # 定義自定義信號，傳遞 ligand_name 和選定的 mode
    affinity_changed = pyqtSignal(str, int, dict)  # 信號屬性，它是 pyqtSignal 類型，會成為每個實例的屬性
    
    def __init__(self, ligand_name, result_data_dict, parent=None):
        super().__init__()
        
        # 獲取 mode 和 affinity 資料
        self.ligand_name = ligand_name
        self.result_data_dict = result_data_dict
        self.mode_list = result_data_dict['mode']
        self.affinity_list = result_data_dict['affinity']
        
        
        # 構建 mode 與 affinity 的字串格式
        self.mode_affinities = [f"{mode}: {affinity:.4f}" for mode, affinity in zip(self.mode_list, self.affinity_list)]
        
        
        # 創建 QComboBox 顯示模式與對應的 affinity 值
        self.combo_box = QComboBox()
        self.combo_box.addItems(self.mode_affinities)

        # 創建左右切換按鈕
        self.left_button = QPushButton("<")
        self.right_button = QPushButton(">")
        self.left_button.setFixedWidth(20)
        self.right_button.setFixedWidth(20)

        # 佈局設置
        layout = QHBoxLayout()
        layout.addWidget(self.left_button)
        layout.addWidget(self.combo_box)
        layout.addWidget(self.right_button)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        # 設置初始值
        self.combo_box.setCurrentIndex(0)
        
        # 信號連接
        self.left_button.clicked.connect(self.move_left)
        self.right_button.clicked.connect(self.move_right)
        
        # 當選擇模式改變時發送信號
        self.combo_box.currentIndexChanged.connect(self.emit_affinity_changed_signal)

    def move_left(self):
        """切換到前一個 mode"""
        current_index = self.combo_box.currentIndex()
        if current_index > 0:
            self.combo_box.setCurrentIndex(current_index - 1)

    def move_right(self):
        """切換到下一個 mode"""
        current_index = self.combo_box.currentIndex()
        if current_index < len(self.mode_affinities) - 1:
            self.combo_box.setCurrentIndex(current_index + 1)
    
    def emit_affinity_changed_signal(self):
        """發送選擇的 mode 改變信號"""
        selected_mode = self.combo_box.currentIndex() + 1  # mode 索引從 1 開始
        self.affinity_changed.emit(self.ligand_name, selected_mode, self.result_data_dict)  # 傳遞 ligand 名稱和模式

    def get_current_affinity(self):
        """返回當前選中的 mode 與 affinity 值"""
        return self.combo_box.currentText()
    
    
    

    
    
    
        