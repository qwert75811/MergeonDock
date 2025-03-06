rgeonDock - Open-Source Molecular Docking & Analysis Tool  

MergeonDock is an open-source tool for molecular docking and interaction analysis, integrating AutoDock, RDKit, PyMOL, and Open Babel to streamline docking workflows.

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

---

## **Platform Compatibility**
MergeonDock is developed primarily for **Windows**.

### **Windows**
- Uses third-party Windows wheels for PyMOL and Open Babel.
- The **installer includes these dependencies**, so users do not need to install them manually.

### **Linux**
- **The current codebase is not directly compatible with Linux.**  
- However, developers can modify it to support Linux.  
- PyMOL and Open Babel can be installed using official Linux packages


## **Dependency Handling in MergeonDock**
1️. API-Based Dependencies (import in Python)
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
All trademarks and software rights belong to their respective owners.