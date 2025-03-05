# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:46:40 2024

# MergeonDock - Open-Source Molecular Docking & Analysis Tool  

Copyright (c) 2024 Xhamrock Studio  

Permission is hereby granted, free of charge, to any person obtaining a copy  
of this software and associated documentation files (the "Software"), to deal  
in the Software without restriction, including without limitation the rights  
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell  
copies of the Software, and to permit persons to whom the Software is  
furnished to do so, subject to the following conditions:  

The above copyright notice and this permission notice shall be included in all  
copies or substantial portions of the Software.  

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE  
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER  
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  
SOFTWARE.  

---

## **Technologies Used**
MergeonDock integrates the following open-source technologies:

| Software         | Version  | License | Source Link |
|-----------------|---------|---------|-------------|
| **AutoDock 4.2** | 4.2     | AutoDock License | [AutoDock Website](http://autodock.scripps.edu/) |
| **AutoDock Vina** | 1.2.3  | Apache 2.0 | [AutoDock Vina GitHub](https://github.com/ccsb-scripps/AutoDock-Vina) |
| **RDKit**       | Latest  | BSD | [RDKit Website](https://www.rdkit.org/) |
| **PyMOL (Open-Source)** | 3.1.0 | GPL | [PyMOL GitHub (official)](https://github.com/schrodinger/pymol-open-source) |
| **Open Babel**  | 3.1.1   | GPL/LGPL | [Open Babel GitHub (official)](https://github.com/openbabel/openbabel) |

### **Windows-Specific Third-Party Wheels**
Since PyMOL and Open Babel do not officially provide Python Wheels for Windows, we use third-party compiled binaries:

| Software         | Windows Wheel Source | Maintainer |
|-----------------|----------------------|------------|
| **PyMOL (Windows Wheel)** | [cgohlke/pymol-open-source-wheels](https://github.com/cgohlke/pymol-open-source-wheels) | Third-Party |
| **Open Babel (Windows Wheel)** | [njzjz/openbabel-wheel](https://github.com/njzjz/openbabel-wheel) | Third-Party |

**Note:**  
- These Windows wheels are **not officially maintained** by PyMOL or Open Babel developers.  
- MergeonDock **does not modify these components**; they are included as dependencies.  

---

## **Platform Compatibility**
MergeonDock is developed primarily for **Windows**.

### **Windows**
- Uses third-party Windows wheels for PyMOL and Open Babel.
- The **installer includes these dependencies**, so users do not need to install them manually.

### **Linux**
- **The current codebase is not directly compatible with Linux.**  
- However, developers can modify it to support Linux.  
- PyMOL and Open Babel can be installed using official Linux packages:
  ```bash
  # Install PyMOL (official)
  sudo apt install pymol

  # Install Open Babel (official)
  sudo apt install openbabel


## **Dependency Handling in MergeonDock**
1️.API-Based Dependencies (import in Python)
 PyMOL & Open Babel: Used via Python API (import pymol / import openbabel) without modification.

2️. External Process Execution (QProcess)
 AutoDock Vina:
 The vina.exe binary is included in the project.
 
 MergeonDock calls it externally using QProcess, with command-line arguments provided by the software.
 AutoDock 4.2:

 AutoDock 4.2 is a Python 2-based tool, which cannot run directly in Python 3.
 To resolve this, MergeonDock includes a minimal custom Python 2 environment, containing only the required autodocktools components.
 MergeonDock calls this environment externally via QProcess to execute docking commands.

AutoDock Integration
MergeonDock includes AutoDock binaries for convenience. These binaries are sourced from:

AutoDock 4.2: AutoDock Official Website
AutoDock Vina: AutoDock Vina GitHub
MergeonDock does not modify AutoDock. The software simply calls these tools externally.

Final Notes
MergeonDock is an independent project and is not affiliated with PyMOL, Open Babel, or AutoDock developers.
This software includes GPL-licensed components (PyMOL, Open Babel).
All trademarks and software rights belong to their respective owners.



@author: Xhamrock Studio
"""

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from PyQt5 import QtWidgets

from pmg_qt.pymol_gl_widget import PyMOLGLWidget

from MergeonDock import gui


from MergeonDock.menu import File_format_convert
from MergeonDock.menu import advance_setting
from MergeonDock.menu import about

from MergeonDock import all_parameters
from MergeonDock import work_directory
from MergeonDock.ligands_upload import ligands_upload
from MergeonDock.receptor_upload import file_upload
from MergeonDock.receptor_upload import rec_prepare_detect
from MergeonDock.dock_setting import gridbox
from MergeonDock.dock_setting import parameters
from MergeonDock import dock
from MergeonDock.dock_analysis import dock_analysis_basic


class MyPyMOLGLWidget(PyMOLGLWidget):
    def __init__(self, parent=None):
        super(MyPyMOLGLWidget, self).__init__(parent)
    
    def dragEnterEvent(self, event):
        event.ignore()  # 不接受拖放事件

    def dropEvent(self, event):
        event.ignore()  # 不处理拖放事件



class Start(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = gui.Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.ui.actionQuit.triggered.connect(self.close_program)
       

        self.pymol_process= MyPyMOLGLWidget(self.ui.widget_pymol)
        layout = QGridLayout(self.ui.widget_pymol)  #重新設置布局
        layout.setContentsMargins(0, 0, 0, 0)  # 设置布局边距为0
        layout.setSpacing(0)  # 设置布局间距为0
        layout.addWidget(self.pymol_process)  # 把PyMOL widget添加到布局中
        self.ui.widget_pymol.setLayout(layout)  # 應用布局
        
        
        
        self.current_tab = "dock"
        self.pymol_session = {"dock tab": "", "analysis tab": ""}

        self.setup()
        
        # 连接 Tab 切换信号到处理函数
        self.ui.tabWidget.currentChanged.connect(self.tab_changed)
        
        
        current_widget = self.ui.stackedWidget.currentWidget()
        if current_widget == self.ui.page_show_all_input:
            
            self.ui.pushButton_setgridbox.clicked.connect(self.gridbox_instance.switch_to_gridbox_tab)
            self.ui.pushButton_setparameter.clicked.connect(self.parameters_instance.switch_to_parameter_tab)
        else:
            self.ui.stackedWidget.setCurrentWidget(self.ui.page_show_all_input)

  
    def tab_changed(self):
        if self.current_tab == "dock":
            dock_session = self.pymol_process.cmd.get_session()
            self.pymol_session["dock tab"] = dock_session
            if self.pymol_session["analysis tab"] != "":
                self.pymol_process.cmd.set_session(self.pymol_session["analysis tab"])
            else:
                self.pymol_process.cmd.reinitialize()
            self.current_tab = "analysis"
        elif self.current_tab == "analysis":
            analysis_session = self.pymol_process.cmd.get_session()
            self.pymol_session["analysis tab"] = analysis_session
            if self.pymol_session["dock tab"] != "":
                self.pymol_process.cmd.set_session(self.pymol_session["dock tab"])
            else:
                self.pymol_process.cmd.reinitialize()
            self.current_tab = "dock"
     
        
     
    def setup(self):    
        self.all_parameters = all_parameters.Parameters_storage()
      
        
        #把UI當作參數傳遞給work_directory.py
        self.work_directory_instance = work_directory.Directory_setup(self.ui, self.pymol_process, self.all_parameters)
        
        #把UI當作參數傳遞給ligands_upload.py
        self.ligands_upload_instance = ligands_upload.Ligands_upload(self.ui, self.pymol_process, self.all_parameters)
        
        #把UI當作參數傳遞給receptor_upload.py
        self.receptor_upload_instance = file_upload.Receptor_upload(self.ui, self.pymol_process, self.all_parameters)
        self.receptor_detection_instance = rec_prepare_detect.Receptor_sequence_detection(self.pymol_process, self.all_parameters, self.receptor_upload_instance)
        
        
        #把UI當作參數傳遞給gridbox.py
        self.gridbox_instance = gridbox.Gridbox_setting(self.ui, self.pymol_process, self.all_parameters)

        #把UI當作參數傳遞給parameters.py
        self.parameters_instance = parameters.Parameter_setting(self.ui, self.pymol_process, self.all_parameters)
        
        #互相傳遞實例
        self.gridbox_instance.set_parameters_instance(self.parameters_instance)
        self.parameters_instance.set_gridbox_instance(self.gridbox_instance)
        
        
        
        
        
        
        #把UI當作參數傳遞給Option_format_convert.py
        self.option_format_convert = File_format_convert.Menu_option_file_convert(self.ui, self.all_parameters)
        
        #Advance Setting
        self.option_advence_setting = advance_setting.Menu_option_advance_setting(self.ui, self.all_parameters)
        
        #about
        self.help_about = about.Menu_help_about(self.ui, self.all_parameters)
        
        #分析頁面實例
        self.analysis_basic_tab = dock_analysis_basic.Analysis_results(self.ui, self.pymol_process, self.all_parameters)
        
        #把UI當作參數傳遞給dock.py
        self.dock_instance = dock.Dock_setting(self.ui, self.all_parameters, self.analysis_basic_tab)
        
        
        print("install path:", self.all_parameters.install_path)
    
    
    def close_program(self):
        self.close()
    
    
        
        
    def run(self):
        self.show()


def main():
    app = QApplication(sys.argv)
    window = Start()
    window.run()
    sys.exit(app.exec_())
  
if __name__ == "__main__":
    main()
    
