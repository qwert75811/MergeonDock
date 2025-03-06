---
title: MergeonDock - Open-Source Molecular Docking & Analysis Tool

---

![MergeonDock_title](https://hackmd.io/_uploads/rJUwhP8oJg.png)
# MergeonDock - Open-Source Molecular Docking & Analysis Tool  
MergeonDock is an open-source tool for molecular docking and interaction analysis, integrating AutoDock, RDKit, PyMOL, and Open Babel to streamline docking workflows.

## **Features**
- **Graphical User Interface (GUI):** Provides an intuitive interface for docking workflow visualization.
- **AutoDock Vina Integration:** Uses AutoDock Vina 1.2.5 as the core docking engine.
- **Receptor & Ligand Processing:** Automatically prepares `.pdb` files for docking.
- **Advanced Settings:** Customize receptor and ligand preparation settings based on AutoDock 4.2 options.
- **Mini Open Babel Converter**: A built-in tool for converting molecular file formats between .mol, .pdb, and .pdbqt using Open Babel, seamlessly integrated into the GUI.
- **Interaction Analysis:** Allows users to analyze docking results with reference ligands.
- **Multi-Format Support:** Supports `.pdb`, `.pdbqt`, and `.sdf` file formats.

## **Installation**
**Download and install from my release (Windows OS only).**  
Linux is currently not supported. For details, please refer to the Platform Compatibility section below.

---

## **Development Status**
- This is an early-stage release, and many features are only placeholders without actual functionality. It is normal to encounter issues while using it.
- As I am not a professional software developer, maintaining this project requires significant effort, and there are also potential infringement risks. While future updates are unlikely, they are not impossible. My main goal is education and promotion, which is why I have decided to make all information and code publicly available, allowing those who are interested to explore and expand upon it freely.

## Quick User Guide
This section provides a step-by-step guide to using MergeonDock.

### Step 1. Choose your folder
Each docking run is assigned its own independent folder. Repeating different docking tasks in the same folder may risk overwriting files. You can go back to Step 1 and click 'New Project' to select a new folder.
![step1](https://hackmd.io/_uploads/rk_Wa58iye.jpg)

### Step 2. Upload your Receptor file and Ligands files
It will automatically prepare your raw files into AutoDock PDBQT files.
![step2](https://hackmd.io/_uploads/SykMaqLjye.jpg)

#### Step 2-1. Advanced Setting (Optional, for customization)
This feature is based on AutoDock 4's preparation functions, providing a parameter selection table. Move your cursor over a parameter and hover for a moment to see its original description.

This table is designed to simplify parameter configuration. However, if errors or crashes occur after using it, you may need to test the settings directly in AutoDockTools. This feature only calls AutoDock externally.

This step is not mandatory and can be skipped for general use.
![step2-1](https://hackmd.io/_uploads/SJHrYoLiyl.jpg)

#### Step 2-2. Receptor detection
If your receptor file is in .pdb format and comes from the PDB bank or is a pre-processed file, MergeonDock provides a simple heteroatom recognition mechanism for you to choose from. If you find that the uploaded file does not behave as expected or does not trigger this mechanism, you may need to manually perform the preprocessing.
![step2-2](https://hackmd.io/_uploads/BkOgvs8i1g.jpg)



#### Step 2-3. Ligands support
The current version supports .pdb, .pdbqt, and .sdf files from the database, but preparation may take some time.
![step2-2](https://hackmd.io/_uploads/BJ2Qa5Lsyx.jpg)

### Step 3. Gridbox Setting
Adjust your gridbox parameters (save for later return).
![step3](https://hackmd.io/_uploads/HJBE6qUjye.jpg)


### Step 4. Basic docking parameter setting
Adjust your docking parameters (save for later return).
![step4](https://hackmd.io/_uploads/SykSTqLske.jpg)


### Step 5. Start docking
Press the docking button to start the docking process (this may take some time depending on your CPU).
![step5](https://hackmd.io/_uploads/r1OH65Ljyg.jpg)


### Step 6. Move to the analysis window to view the results
![step6](https://hackmd.io/_uploads/ByxIa9IsJx.jpg)


### Step 7. Interaction analysis
![step7](https://hackmd.io/_uploads/HyII6qUoJl.jpg)



---

## **Technologies Used**
MergeonDock integrates the following open-source technologies:

| Software         | Version  | License | Source Link |
|-----------------|---------|---------|-------------|
| **AutoDock 4.2** | 4.2     | AutoDock License | [AutoDock Website](http://autodock.scripps.edu/) |
| **AutoDock Vina** | 1.2.5  | Apache 2.0 | [AutoDock Vina GitHub](https://github.com/ccsb-scripps/AutoDock-Vina) |
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

## **Dependency Handling in MergeonDock**

### **1️⃣ API-Based Dependencies (`import` in Python)**
✅ **PyMOL** & **Open Babel**:  
- Used via Python API (`import pymol` / `import openbabel`) without modification.  
- These dependencies are installed via Python Wheels.

### **2️⃣ External Process Execution (QProcess)**
✅ **AutoDock Vina**:  
- The `vina.exe` binary is **included in the project**.  
- MergeonDock calls it **externally using QProcess**, with command-line arguments provided by the software.  

✅ **AutoDock 4.2**:  
- AutoDock 4.2 is a **Python 2-based tool**, which **cannot run directly in Python 3**.  
- To resolve this, MergeonDock **includes a minimal custom Python 2 environment**, containing only the required `autodocktools` components.  
- MergeonDock calls this environment **externally via QProcess** to execute docking commands.

These binaries are sourced from:
- **AutoDock 4.2**: [AutoDock Official Website](http://autodock.scripps.edu/)
- **AutoDock Vina**: [AutoDock Vina GitHub](https://github.com/ccsb-scripps/AutoDock-Vina)


## **Platform Compatibility**
MergeonDock is developed primarily for **Windows**.

### **Windows**
- Fully supported. The installer includes all necessary dependencies.

### **Linux**
⚠ **Linux compatibility is currently untested.**  
MergeonDock is designed for Windows, and running it on Linux may require significant modifications.  
Potential issues include:
- **AutoDock Vina & AutoDock 4.2 execution** (command differences & Python 2 environment)
- **PyQt5 GUI behavior** (file dialogs, QProcess handling)
- **PyMOL & Open Babel dependencies** (must use official Linux versions)

Developers interested in porting MergeonDock to Linux may need to manually adjust these components.

## **License**
MergeonDock is released under the **MIT License**.

MIT License

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


