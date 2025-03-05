# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:40:08 2024

@author: Xhamrock Studio
"""


import os, re
from PyQt5.QtWidgets import QDialog
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QTableWidgetItem, QRadioButton, QButtonGroup, QWidget, QVBoxLayout
from PyQt5.QtCore import Qt
import subprocess

from MergeonDock.dock_setting import choose_ref_ligand_ui
from MergeonDock.dock_setting import pymol_gridbox




class Gridbox_setting():
    def __init__(self, ui, pymol_process, all_parameters):
        self.ui = ui
        self.pymol_process = pymol_process  # 將 pymol_process 參數保存為類的屬性
        self.all_parameters = all_parameters

        self.ui.pushButton_use_ref_ligands.clicked.connect(self.ref_ligand_for_grid)

        # Connect sliders and spinboxes
        self.ui.horizontalSlider_centerX.valueChanged.connect(self.update_spinbox_centerX)
        self.ui.doubleSpinBox_centerX.valueChanged.connect(self.update_slider_centerX)
        
        self.ui.horizontalSlider_centerY.valueChanged.connect(self.update_spinbox_centerY)
        self.ui.doubleSpinBox_centerY.valueChanged.connect(self.update_slider_centerY)
        
        self.ui.horizontalSlider_centerZ.valueChanged.connect(self.update_spinbox_centerZ)
        self.ui.doubleSpinBox_centerZ.valueChanged.connect(self.update_slider_centerZ)

        self.ui.horizontalSlider_sizeX.valueChanged.connect(self.update_spinbox_sizeX)
        self.ui.spinBox_sizeX.valueChanged.connect(self.update_slider_sizeX)
        
        self.ui.horizontalSlider_sizeY.valueChanged.connect(self.update_spinbox_sizeY)
        self.ui.spinBox_sizeY.valueChanged.connect(self.update_slider_sizeY)
        
        self.ui.horizontalSlider_sizeZ.valueChanged.connect(self.update_spinbox_sizeZ)
        self.ui.spinBox_sizeZ.valueChanged.connect(self.update_slider_sizeZ)
        
        self.ui.doubleSpinBox_spacing.valueChanged.connect(self.update_space)
    
    def set_parameters_instance(self, parameters_instance):
        self.parameters_instance = parameters_instance       
        
    
    def update_spinbox_centerX(self, value):
        self.ui.doubleSpinBox_centerX.setValue(value / 1000.0)
        self.gridbox.update_center(self.ui.doubleSpinBox_centerX.value(),
                                   self.ui.doubleSpinBox_centerY.value(),
                                   self.ui.doubleSpinBox_centerZ.value())
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
        
    def update_slider_centerX(self, value):
        self.ui.horizontalSlider_centerX.setValue(int(value * 1000))
        
    def update_spinbox_centerY(self, value):
        self.ui.doubleSpinBox_centerY.setValue(value / 1000.0)
        self.gridbox.update_center(self.ui.doubleSpinBox_centerX.value(),
                                   self.ui.doubleSpinBox_centerY.value(),
                                   self.ui.doubleSpinBox_centerZ.value())
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
        
    def update_slider_centerY(self, value):
        self.ui.horizontalSlider_centerY.setValue(int(value * 1000))
        
    def update_spinbox_centerZ(self, value):
        self.ui.doubleSpinBox_centerZ.setValue(value / 1000.0)
        self.gridbox.update_center(self.ui.doubleSpinBox_centerX.value(),
                                   self.ui.doubleSpinBox_centerY.value(),
                                   self.ui.doubleSpinBox_centerZ.value())
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
        
    def update_slider_centerZ(self, value):
        self.ui.horizontalSlider_centerZ.setValue(int(value * 1000))
        
    def update_spinbox_sizeX(self, value):
        self.ui.spinBox_sizeX.setValue(value)
        self.gridbox.update_size(self.ui.spinBox_sizeX.value(),
                                 self.ui.spinBox_sizeY.value(),
                                 self.ui.spinBox_sizeZ.value())
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
        
    def update_slider_sizeX(self, value):
        self.ui.horizontalSlider_sizeX.setValue(value)
        
    def update_spinbox_sizeY(self, value):
        self.ui.spinBox_sizeY.setValue(value)
        self.gridbox.update_size(self.ui.spinBox_sizeX.value(),
                                 self.ui.spinBox_sizeY.value(),
                                 self.ui.spinBox_sizeZ.value())
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
        
    def update_slider_sizeY(self, value):
        self.ui.horizontalSlider_sizeY.setValue(value)
        
    def update_spinbox_sizeZ(self, value):
        self.ui.spinBox_sizeZ.setValue(value)
        self.gridbox.update_size(self.ui.spinBox_sizeX.value(),
                                 self.ui.spinBox_sizeY.value(),
                                 self.ui.spinBox_sizeZ.value())
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
        
    def update_slider_sizeZ(self, value):
        self.ui.horizontalSlider_sizeZ.setValue(value)
    
    def update_space(self, value):
        self.gridbox.update_space(value)
        self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
   
    
    
    
    def switch_to_gridbox_tab(self):
        # 解绑之前的连接
        self.ui.pushButton_setgridbox.clicked.disconnect()
        self.ui.pushButton_setparameter.clicked.disconnect()
    
        # 绑定 gridbox 特定的功能
        self.ui.pushButton_setgridbox.clicked.connect(self.grid_button_save_button)
        self.ui.pushButton_setparameter.clicked.connect(self.parameter_button_cancel_button)
    
        # 更新按钮文本
        self.ui.pushButton_setgridbox.setText("Save")
        self.ui.pushButton_setparameter.setText("Cancel")
        
        #切換頁面
        self.ui.stackedWidget.setCurrentWidget(self.ui.page_gridbox_setting)
       
        
        current_widget = self.ui.stackedWidget.currentWidget()
        if current_widget == self.ui.page_gridbox_setting:
            self.pymol_process.cmd.zoom(self.all_parameters.output_prepared_receptor_name) #畫面鎖在receptor
            # 根据条件创建 gridbox 实例
            if [self.ui.doubleSpinBox_centerX.value(), self.ui.doubleSpinBox_centerY.value(), self.ui.doubleSpinBox_centerZ.value()] == [0,0,0]:
                center = self.pymol_process.cmd.centerofmass(self.all_parameters.output_prepared_receptor_name)                  #抓取receptor的中心值
                size = [self.ui.spinBox_sizeX.value(), self.ui.spinBox_sizeY.value(), self.ui.spinBox_sizeZ.value()]    
                space = round(self.ui.doubleSpinBox_spacing.value(), 3)
                self.gridbox = pymol_gridbox.PyMOLGridBox(center=center, size=size, space=space)
                
                #套用receptor中心值當預設
                self.ui.doubleSpinBox_centerX.setValue(center[0])
                self.ui.doubleSpinBox_centerY.setValue(center[1])
                self.ui.doubleSpinBox_centerZ.setValue(center[2])
                self.all_parameters.gridcenter_X = self.ui.doubleSpinBox_centerX.value()
                self.all_parameters.gridcenter_Y = self.ui.doubleSpinBox_centerY.value()
                self.all_parameters.gridcenter_Z = self.ui.doubleSpinBox_centerZ.value()

            else:
                center = [self.ui.doubleSpinBox_centerX.value(), self.ui.doubleSpinBox_centerY.value(), self.ui.doubleSpinBox_centerZ.value()]  
                size = [self.ui.spinBox_sizeX.value(), self.ui.spinBox_sizeY.value(), self.ui.spinBox_sizeZ.value()] 
                space = round(self.ui.doubleSpinBox_spacing.value(), 3)
                self.gridbox = pymol_gridbox.PyMOLGridBox(center=center, size=size, space=space)

                self.pymol_process.cmd.load_cgo(self.gridbox.draw_colored_box(), 'gridbox')
    
    
    
    
    def grid_button_save_button(self):
        self.pymol_process.cmd.delete("gridbox")
        
        centerX_value = self.ui.doubleSpinBox_centerX.value()
        centerY_value = self.ui.doubleSpinBox_centerY.value()
        centerZ_value = self.ui.doubleSpinBox_centerZ.value()
        sizeX_value = self.ui.spinBox_sizeX.value()
        sizeY_value = self.ui.spinBox_sizeY.value()
        sizeZ_value = self.ui.spinBox_sizeZ.value()
        spacing_value = round(self.ui.doubleSpinBox_spacing.value(), 3)
      
        self.all_parameters.gridcenter_X = centerX_value
        self.all_parameters.gridcenter_Y = centerY_value
        self.all_parameters.gridcenter_Z = centerZ_value
        self.all_parameters.gridsize_X = sizeX_value
        self.all_parameters.gridsize_Y = sizeY_value
        self.all_parameters.gridsize_Z = sizeZ_value
        self.all_parameters.spacing_value = spacing_value
        
        
        self.ui.stackedWidget.setCurrentWidget(self.ui.page_show_all_input)
        self.ui.pushButton_setgridbox.setText("Grid box")
        self.ui.pushButton_setparameter.setText("Parameters")
        
        self.ui.pushButton_setgridbox.clicked.disconnect()
        self.ui.pushButton_setparameter.clicked.disconnect()
        
        self.ui.pushButton_setgridbox.clicked.connect(self.switch_to_gridbox_tab)
        self.ui.pushButton_setparameter.clicked.connect(self.parameters_instance.switch_to_parameter_tab)
        

    def parameter_button_cancel_button(self):
        current_widget = self.ui.stackedWidget.currentWidget()
        if current_widget == self.ui.page_gridbox_setting:
            self.ui.doubleSpinBox_centerX.setValue(self.all_parameters.gridcenter_X)
            self.ui.doubleSpinBox_centerY.setValue(self.all_parameters.gridcenter_Y)
            self.ui.doubleSpinBox_centerZ.setValue(self.all_parameters.gridcenter_Z)
            self.ui.spinBox_sizeX.setValue(self.all_parameters.gridsize_X)
            self.ui.spinBox_sizeY.setValue(self.all_parameters.gridsize_Y)
            self.ui.spinBox_sizeZ.setValue(self.all_parameters.gridsize_Z)
            self.ui.doubleSpinBox_spacing.setValue(self.all_parameters.spacing_value)
            self.pymol_process.cmd.delete("gridbox")
            self.ui.stackedWidget.setCurrentWidget(self.ui.page_show_all_input)
            self.ui.pushButton_setgridbox.setText("Grid box")
            self.ui.pushButton_setparameter.setText("Parameters")
        
            self.ui.pushButton_setgridbox.clicked.disconnect()
            self.ui.pushButton_setparameter.clicked.disconnect()
            self.ui.pushButton_setgridbox.clicked.connect(self.switch_to_gridbox_tab)
            self.ui.pushButton_setparameter.clicked.connect(self.parameters_instance.switch_to_parameter_tab)
            
       
    
    def ref_ligand_for_grid(self):
        if len(self.all_parameters.ref_prepared_ligands_path) > 1:
            self.pick_window_ui = choose_ref_ligand_ui.Ui_Dialog_pick_ref_ligand()
            self.pick_window = QDialog()
            self.pick_window_ui.setupUi(self.pick_window)
            self.pick_window.show()
            
            self.pick_window_ui.tableWidget_reflig_list.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            
            self.pick_window_ui.pushButton_open.clicked.connect(self.open_button)
            self.pick_window_ui.pushButton_confirm.clicked.connect(self.confirm_button)
            self.pick_window_ui.pushButton_cancel.clicked.connect(self.cancel_button)
            
            self.pick_window_ui.tableWidget_reflig_list.setRowCount(len(self.all_parameters.ref_prepared_ligands_path))
            for index in range(int(len(self.all_parameters.ref_prepared_ligands_name))):
                self.pick_window_ui.tableWidget_reflig_list.setItem(index, 0, QTableWidgetItem(self.all_parameters.ref_prepared_ligands_name[index]))
            
            
            self.radio_buttongroup = QButtonGroup()
            for row in range(int(len(self.all_parameters.ref_prepared_ligands_name))):
                temp_radiobutton = QRadioButton()
                self.radio_buttongroup.addButton(temp_radiobutton)
                temp_widget = QWidget()
                temp_layout = QVBoxLayout()
                temp_layout.addWidget(temp_radiobutton)
                temp_layout.setAlignment(temp_radiobutton, Qt.AlignCenter)
                temp_layout.setContentsMargins(0, 0, 0, 0)
                temp_widget.setLayout(temp_layout)
                
                self.pick_window_ui.tableWidget_reflig_list.setCellWidget(row, 1, temp_widget)
                            
     
        elif len(self.all_parameters.ref_prepared_ligands_path) == 1:
            try:
                self.output_gpf_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.output_prepared_receptor_name}_ref_gpf"))
                self.all_parameters.ref_ligand_picked_path = self.all_parameters.ref_prepared_ligands_path[0]
                subprocess.run(f'{self.all_parameters.autodock4_run_prepare_gpf} -l "{self.all_parameters.ref_ligand_picked_path}" -r "{self.all_parameters.output_prepared_receptor_path}" -o "{self.output_gpf_path}" -y')
                if os.path.exists(self.output_gpf_path):
                    self.obtain_value_from_auto_gpf()
                else:
                    print("auto gpf file not found")
            except:
                print("fail")
        else:
            ref_ligand_path = QtWidgets.QFileDialog.getOpenFileName(None, "Choose Ref. ligand", "", "pdbqt files (*.pdbqt)" )
            if os.path.isfile(ref_ligand_path[0]):
                ref_ligand_upload_path = os.path.normpath(ref_ligand_path[0])
                
                self.output_gpf_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.output_prepared_receptor_name}_ref_gpf"))
                subprocess.run(f'{self.all_parameters.autodock4_run_prepare_gpf} -l "{ref_ligand_upload_path}" -r "{self.all_parameters.output_prepared_receptor_path}" -o "{self.output_gpf_path}" -y')
                
                if os.path.exists(self.output_gpf_path):
                    self.obtain_value_from_auto_gpf()
                else:
                    print("auto gpf file not found")
            else:
                print("No ref was picked")
    
    
    def open_button(self):
        ref_ligand_path = QtWidgets.QFileDialog.getOpenFileName(None, "Choose Ref. ligand", "", "pdbqt files (*.pdbqt)" )
        self.all_parameters.ref_ligand_picked_path = os.path.normpath(ref_ligand_path[0])
        
        self.output_gpf_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.output_prepared_receptor_name}_ref_gpf"))
        subprocess.run(f'{self.all_parameters.autodock4_run_prepare_gpf} -l "{self.all_parameters.ref_ligand_picked_path}" -r "{self.all_parameters.output_prepared_receptor_path}" -o "{self.output_gpf_path}" -y')
        
        if os.path.exists(self.output_gpf_path):
            self.obtain_value_from_auto_gpf()
        else:
            print("auto gpf file not found")
        
        self.pick_window.close()
            
    
    def confirm_button(self):
        selected_button = self.radio_buttongroup.checkedButton()
        if selected_button:
            selected_row = self.pick_window_ui.tableWidget_reflig_list.indexAt(selected_button.pos()).row()
            self.all_parameters.ref_ligand_picked_path = self.all_parameters.ref_prepared_ligands_path[selected_row]

            self.output_gpf_path = os.path.normpath(os.path.join(self.all_parameters.work_directory, f"{self.all_parameters.output_prepared_receptor_name}_ref_gpf"))
            subprocess.run(f'{self.all_parameters.autodock4_run_prepare_gpf} -l "{self.all_parameters.ref_ligand_picked_path}" -r "{self.all_parameters.output_prepared_receptor_path}" -o "{self.output_gpf_path}" -y')

            if os.path.exists(self.output_gpf_path):
                self.obtain_value_from_auto_gpf()
            else:
                print("auto gpf file not found")

        self.pick_window.close()
              
       
    def cancel_button(self):
        self.pick_window.close()
            
        
    def obtain_value_from_auto_gpf(self):
        self.all_parameters.gpf_file_path = self.output_gpf_path
        with open(self.all_parameters.gpf_file_path, "r") as gpf_file:
            self.full_content = gpf_file.read()
            self.full_line_content = self.full_content.splitlines()
        
        pattern_npt = re.compile(r'npts\s+(\d+)\s+(\d+)\s+(\d+)\s+')
        for line in self.full_line_content:
            match = pattern_npt.match(line)
            if match:
                npt_x = match.group(1) 
                npt_y = match.group(2)
                npt_z = match.group(3)
                break
                
        self.ui.spinBox_sizeX.setValue(int(npt_x))
        self.ui.spinBox_sizeY.setValue(int(npt_y))
        self.ui.spinBox_sizeZ.setValue(int(npt_z))
        self.ui.horizontalSlider_sizeX.setValue(int(npt_x))
        self.ui.horizontalSlider_sizeY.setValue(int(npt_y))
        self.ui.horizontalSlider_sizeY.setValue(int(npt_z))
         
        pattern_gridcenter = re.compile(r'gridcenter\s+([+-]?\d*\.\d+)\s+([+-]?\d*\.\d+)\s+([+-]?\d*\.\d+)\s+')
        for line in self.full_line_content:
            match = pattern_gridcenter.match(line)
            if match:
                center_x = match.group(1) 
                center_y = match.group(2)
                center_z = match.group(3)
                break
                
        self.ui.doubleSpinBox_centerX.setValue(float(center_x))
        self.ui.doubleSpinBox_centerY.setValue(float(center_y))
        self.ui.doubleSpinBox_centerZ.setValue(float(center_z))
        self.ui.horizontalSlider_centerX.setValue(int(float(center_x) * 1000))
        self.ui.horizontalSlider_centerY.setValue(int(float(center_y) * 1000))
        self.ui.horizontalSlider_centerZ.setValue(int(float(center_z) * 1000))
        
                
       
            
    
    
        