# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'advance_setting_ui.ui'
#
# Created by: PyQt5 UI code generator 5.15.11
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog_advance_setting(object):
    def setupUi(self, Dialog_advance_setting):
        Dialog_advance_setting.setObjectName("Dialog_advance_setting")
        Dialog_advance_setting.resize(686, 680)
        self.gridLayout_8 = QtWidgets.QGridLayout(Dialog_advance_setting)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.groupBox_ligand_preparation = QtWidgets.QGroupBox(Dialog_advance_setting)
        self.groupBox_ligand_preparation.setStyleSheet("font: 75 12pt \"Arial\";")
        self.groupBox_ligand_preparation.setObjectName("groupBox_ligand_preparation")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.groupBox_ligand_preparation)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.radioButton_lig_ad4 = QtWidgets.QRadioButton(self.groupBox_ligand_preparation)
        self.radioButton_lig_ad4.setChecked(True)
        self.radioButton_lig_ad4.setObjectName("radioButton_lig_ad4")
        self.gridLayout_7.addWidget(self.radioButton_lig_ad4, 0, 0, 1, 1)
        self.radioButton_lig_meeko = QtWidgets.QRadioButton(self.groupBox_ligand_preparation)
        self.radioButton_lig_meeko.setEnabled(False)
        self.radioButton_lig_meeko.setObjectName("radioButton_lig_meeko")
        self.gridLayout_7.addWidget(self.radioButton_lig_meeko, 0, 1, 1, 1)
        self.stackedWidget_lig_ad4_opt_parameters = QtWidgets.QStackedWidget(self.groupBox_ligand_preparation)
        self.stackedWidget_lig_ad4_opt_parameters.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.stackedWidget_lig_ad4_opt_parameters.setObjectName("stackedWidget_lig_ad4_opt_parameters")
        self.page_lig_ad4 = QtWidgets.QWidget()
        self.page_lig_ad4.setObjectName("page_lig_ad4")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.page_lig_ad4)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.checkBox_lig_ad4_opt_parameters = QtWidgets.QCheckBox(self.page_lig_ad4)
        self.checkBox_lig_ad4_opt_parameters.setObjectName("checkBox_lig_ad4_opt_parameters")
        self.gridLayout_9.addWidget(self.checkBox_lig_ad4_opt_parameters, 0, 0, 1, 2)
        self.splitter_lig_ad4_A = QtWidgets.QSplitter(self.page_lig_ad4)
        self.splitter_lig_ad4_A.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_lig_ad4_A.setObjectName("splitter_lig_ad4_A")
        self.checkBox_lig_ad4_A = QtWidgets.QCheckBox(self.splitter_lig_ad4_A)
        self.checkBox_lig_ad4_A.setEnabled(False)
        self.checkBox_lig_ad4_A.setObjectName("checkBox_lig_ad4_A")
        self.comboBox_lig_ad4_A_setting = QtWidgets.QComboBox(self.splitter_lig_ad4_A)
        self.comboBox_lig_ad4_A_setting.setEnabled(False)
        self.comboBox_lig_ad4_A_setting.setEditable(False)
        self.comboBox_lig_ad4_A_setting.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToContents)
        self.comboBox_lig_ad4_A_setting.setObjectName("comboBox_lig_ad4_A_setting")
        self.comboBox_lig_ad4_A_setting.addItem("")
        self.comboBox_lig_ad4_A_setting.addItem("")
        self.comboBox_lig_ad4_A_setting.addItem("")
        self.comboBox_lig_ad4_A_setting.addItem("")
        self.gridLayout_9.addWidget(self.splitter_lig_ad4_A, 1, 0, 1, 3)
        self.splitter_lig_ad4_R = QtWidgets.QSplitter(self.page_lig_ad4)
        self.splitter_lig_ad4_R.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_lig_ad4_R.setObjectName("splitter_lig_ad4_R")
        self.checkBox_lig_ad4_R = QtWidgets.QCheckBox(self.splitter_lig_ad4_R)
        self.checkBox_lig_ad4_R.setEnabled(False)
        self.checkBox_lig_ad4_R.setMaximumSize(QtCore.QSize(40, 16777215))
        self.checkBox_lig_ad4_R.setObjectName("checkBox_lig_ad4_R")
        self.lineEdit_lig_ad4_R = QtWidgets.QLineEdit(self.splitter_lig_ad4_R)
        self.lineEdit_lig_ad4_R.setEnabled(False)
        self.lineEdit_lig_ad4_R.setMaximumSize(QtCore.QSize(60, 16777215))
        self.lineEdit_lig_ad4_R.setObjectName("lineEdit_lig_ad4_R")
        self.gridLayout_9.addWidget(self.splitter_lig_ad4_R, 3, 0, 1, 1)
        self.splitter_lig_ad4_I = QtWidgets.QSplitter(self.page_lig_ad4)
        self.splitter_lig_ad4_I.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_lig_ad4_I.setObjectName("splitter_lig_ad4_I")
        self.checkBox_lig_ad4_I = QtWidgets.QCheckBox(self.splitter_lig_ad4_I)
        self.checkBox_lig_ad4_I.setEnabled(False)
        self.checkBox_lig_ad4_I.setMaximumSize(QtCore.QSize(40, 16777215))
        self.checkBox_lig_ad4_I.setObjectName("checkBox_lig_ad4_I")
        self.lineEdit_lig_ad4_I = QtWidgets.QLineEdit(self.splitter_lig_ad4_I)
        self.lineEdit_lig_ad4_I.setEnabled(False)
        self.lineEdit_lig_ad4_I.setMaximumSize(QtCore.QSize(200, 16777215))
        self.lineEdit_lig_ad4_I.setObjectName("lineEdit_lig_ad4_I")
        self.gridLayout_9.addWidget(self.splitter_lig_ad4_I, 3, 1, 1, 2)
        self.splitter_lig_ad4_p = QtWidgets.QSplitter(self.page_lig_ad4)
        self.splitter_lig_ad4_p.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_lig_ad4_p.setObjectName("splitter_lig_ad4_p")
        self.checkBox_lig_ad4_p = QtWidgets.QCheckBox(self.splitter_lig_ad4_p)
        self.checkBox_lig_ad4_p.setEnabled(False)
        self.checkBox_lig_ad4_p.setMaximumSize(QtCore.QSize(40, 16777215))
        self.checkBox_lig_ad4_p.setObjectName("checkBox_lig_ad4_p")
        self.lineEdit_lig_ad4_p = QtWidgets.QLineEdit(self.splitter_lig_ad4_p)
        self.lineEdit_lig_ad4_p.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_lig_ad4_p.sizePolicy().hasHeightForWidth())
        self.lineEdit_lig_ad4_p.setSizePolicy(sizePolicy)
        self.lineEdit_lig_ad4_p.setDragEnabled(False)
        self.lineEdit_lig_ad4_p.setObjectName("lineEdit_lig_ad4_p")
        self.gridLayout_9.addWidget(self.splitter_lig_ad4_p, 4, 0, 1, 2)
        self.splitter_lig_ad4_d = QtWidgets.QSplitter(self.page_lig_ad4)
        self.splitter_lig_ad4_d.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_lig_ad4_d.setObjectName("splitter_lig_ad4_d")
        self.checkBox_lig_ad4_d = QtWidgets.QCheckBox(self.splitter_lig_ad4_d)
        self.checkBox_lig_ad4_d.setEnabled(False)
        self.checkBox_lig_ad4_d.setMaximumSize(QtCore.QSize(40, 16777215))
        self.checkBox_lig_ad4_d.setObjectName("checkBox_lig_ad4_d")
        self.lineEdit_lig_ad4_d = QtWidgets.QLineEdit(self.splitter_lig_ad4_d)
        self.lineEdit_lig_ad4_d.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_lig_ad4_d.sizePolicy().hasHeightForWidth())
        self.lineEdit_lig_ad4_d.setSizePolicy(sizePolicy)
        self.lineEdit_lig_ad4_d.setClearButtonEnabled(False)
        self.lineEdit_lig_ad4_d.setObjectName("lineEdit_lig_ad4_d")
        self.gridLayout_9.addWidget(self.splitter_lig_ad4_d, 4, 2, 1, 1)
        self.groupBox_lig_ad4_U_setting = QtWidgets.QGroupBox(self.page_lig_ad4)
        self.groupBox_lig_ad4_U_setting.setEnabled(False)
        self.groupBox_lig_ad4_U_setting.setFlat(False)
        self.groupBox_lig_ad4_U_setting.setCheckable(True)
        self.groupBox_lig_ad4_U_setting.setChecked(False)
        self.groupBox_lig_ad4_U_setting.setObjectName("groupBox_lig_ad4_U_setting")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_lig_ad4_U_setting)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.checkBox_lig_ad4_U_nphs = QtWidgets.QCheckBox(self.groupBox_lig_ad4_U_setting)
        self.checkBox_lig_ad4_U_nphs.setChecked(True)
        self.checkBox_lig_ad4_U_nphs.setObjectName("checkBox_lig_ad4_U_nphs")
        self.gridLayout_5.addWidget(self.checkBox_lig_ad4_U_nphs, 0, 0, 1, 1)
        self.checkBox_lig_ad4_U_lps = QtWidgets.QCheckBox(self.groupBox_lig_ad4_U_setting)
        self.checkBox_lig_ad4_U_lps.setChecked(True)
        self.checkBox_lig_ad4_U_lps.setObjectName("checkBox_lig_ad4_U_lps")
        self.gridLayout_5.addWidget(self.checkBox_lig_ad4_U_lps, 0, 1, 1, 1)
        self.gridLayout_9.addWidget(self.groupBox_lig_ad4_U_setting, 5, 0, 1, 2)
        self.groupBox_lig_ad4_B_setting = QtWidgets.QGroupBox(self.page_lig_ad4)
        self.groupBox_lig_ad4_B_setting.setEnabled(False)
        self.groupBox_lig_ad4_B_setting.setCheckable(True)
        self.groupBox_lig_ad4_B_setting.setChecked(False)
        self.groupBox_lig_ad4_B_setting.setObjectName("groupBox_lig_ad4_B_setting")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.groupBox_lig_ad4_B_setting)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.checkBox_lig_ad4_B_backbone = QtWidgets.QCheckBox(self.groupBox_lig_ad4_B_setting)
        self.checkBox_lig_ad4_B_backbone.setChecked(True)
        self.checkBox_lig_ad4_B_backbone.setObjectName("checkBox_lig_ad4_B_backbone")
        self.gridLayout_6.addWidget(self.checkBox_lig_ad4_B_backbone, 0, 0, 1, 1)
        self.checkBox_lig_ad4_B_amide = QtWidgets.QCheckBox(self.groupBox_lig_ad4_B_setting)
        self.checkBox_lig_ad4_B_amide.setObjectName("checkBox_lig_ad4_B_amide")
        self.gridLayout_6.addWidget(self.checkBox_lig_ad4_B_amide, 0, 1, 1, 1)
        self.checkBox_lig_ad4_B_guanidinium = QtWidgets.QCheckBox(self.groupBox_lig_ad4_B_setting)
        self.checkBox_lig_ad4_B_guanidinium.setObjectName("checkBox_lig_ad4_B_guanidinium")
        self.gridLayout_6.addWidget(self.checkBox_lig_ad4_B_guanidinium, 0, 2, 1, 1)
        self.gridLayout_9.addWidget(self.groupBox_lig_ad4_B_setting, 5, 2, 1, 1)
        self.splitter_lig_ad4_single = QtWidgets.QSplitter(self.page_lig_ad4)
        self.splitter_lig_ad4_single.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_lig_ad4_single.setObjectName("splitter_lig_ad4_single")
        self.checkBox_lig_ad4_C = QtWidgets.QCheckBox(self.splitter_lig_ad4_single)
        self.checkBox_lig_ad4_C.setEnabled(False)
        self.checkBox_lig_ad4_C.setObjectName("checkBox_lig_ad4_C")
        self.checkBox_lig_ad4_Z = QtWidgets.QCheckBox(self.splitter_lig_ad4_single)
        self.checkBox_lig_ad4_Z.setEnabled(False)
        self.checkBox_lig_ad4_Z.setObjectName("checkBox_lig_ad4_Z")
        self.checkBox_lig_ad4_g = QtWidgets.QCheckBox(self.splitter_lig_ad4_single)
        self.checkBox_lig_ad4_g.setEnabled(False)
        self.checkBox_lig_ad4_g.setObjectName("checkBox_lig_ad4_g")
        self.checkBox_lig_ad4_s = QtWidgets.QCheckBox(self.splitter_lig_ad4_single)
        self.checkBox_lig_ad4_s.setEnabled(False)
        self.checkBox_lig_ad4_s.setObjectName("checkBox_lig_ad4_s")
        self.checkBox_lig_ad4_w = QtWidgets.QCheckBox(self.splitter_lig_ad4_single)
        self.checkBox_lig_ad4_w.setEnabled(False)
        self.checkBox_lig_ad4_w.setObjectName("checkBox_lig_ad4_w")
        self.checkBox_lig_ad4_F = QtWidgets.QCheckBox(self.splitter_lig_ad4_single)
        self.checkBox_lig_ad4_F.setEnabled(False)
        self.checkBox_lig_ad4_F.setObjectName("checkBox_lig_ad4_F")
        self.gridLayout_9.addWidget(self.splitter_lig_ad4_single, 2, 0, 1, 2)
        self.stackedWidget_lig_ad4_opt_parameters.addWidget(self.page_lig_ad4)
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.checkBox_meeko_lig_opt_parameters = QtWidgets.QCheckBox(self.page_2)
        self.checkBox_meeko_lig_opt_parameters.setGeometry(QtCore.QRect(30, 20, 191, 16))
        self.checkBox_meeko_lig_opt_parameters.setObjectName("checkBox_meeko_lig_opt_parameters")
        self.stackedWidget_lig_ad4_opt_parameters.addWidget(self.page_2)
        self.gridLayout_7.addWidget(self.stackedWidget_lig_ad4_opt_parameters, 1, 0, 1, 2)
        self.gridLayout_8.addWidget(self.groupBox_ligand_preparation, 1, 0, 1, 1)
        self.groupBox_reduce_onoff = QtWidgets.QGroupBox(Dialog_advance_setting)
        self.groupBox_reduce_onoff.setStyleSheet("font: 75 12pt \"Arial\";")
        self.groupBox_reduce_onoff.setCheckable(True)
        self.groupBox_reduce_onoff.setChecked(False)
        self.groupBox_reduce_onoff.setObjectName("groupBox_reduce_onoff")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox_reduce_onoff)
        self.gridLayout.setObjectName("gridLayout")
        self.radioButton_FILP = QtWidgets.QRadioButton(self.groupBox_reduce_onoff)
        self.radioButton_FILP.setStatusTip("")
        self.radioButton_FILP.setWhatsThis("")
        self.radioButton_FILP.setCheckable(True)
        self.radioButton_FILP.setObjectName("radioButton_FILP")
        self.gridLayout.addWidget(self.radioButton_FILP, 0, 0, 1, 1)
        self.radioButton_NOFLIP = QtWidgets.QRadioButton(self.groupBox_reduce_onoff)
        self.radioButton_NOFLIP.setCheckable(True)
        self.radioButton_NOFLIP.setObjectName("radioButton_NOFLIP")
        self.gridLayout.addWidget(self.radioButton_NOFLIP, 0, 1, 1, 1)
        self.radioButton_Trim = QtWidgets.QRadioButton(self.groupBox_reduce_onoff)
        self.radioButton_Trim.setCheckable(True)
        self.radioButton_Trim.setChecked(True)
        self.radioButton_Trim.setObjectName("radioButton_Trim")
        self.gridLayout.addWidget(self.radioButton_Trim, 0, 2, 1, 1)
        self.radioButton_BUILD = QtWidgets.QRadioButton(self.groupBox_reduce_onoff)
        self.radioButton_BUILD.setCheckable(True)
        self.radioButton_BUILD.setObjectName("radioButton_BUILD")
        self.gridLayout.addWidget(self.radioButton_BUILD, 0, 3, 1, 1)
        self.gridLayout_8.addWidget(self.groupBox_reduce_onoff, 2, 0, 1, 1)
        self.splitter_save_cancel = QtWidgets.QSplitter(Dialog_advance_setting)
        self.splitter_save_cancel.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_save_cancel.setObjectName("splitter_save_cancel")
        self.pushButton_Save = QtWidgets.QPushButton(self.splitter_save_cancel)
        self.pushButton_Save.setMinimumSize(QtCore.QSize(0, 30))
        self.pushButton_Save.setStyleSheet("font: 75 12pt \"Arial\";")
        self.pushButton_Save.setObjectName("pushButton_Save")
        self.pushButton_Cancel = QtWidgets.QPushButton(self.splitter_save_cancel)
        self.pushButton_Cancel.setMinimumSize(QtCore.QSize(0, 30))
        self.pushButton_Cancel.setStyleSheet("font: 75 12pt \"Arial\";")
        self.pushButton_Cancel.setObjectName("pushButton_Cancel")
        self.gridLayout_8.addWidget(self.splitter_save_cancel, 3, 0, 1, 1)
        self.groupBox_receptor_preparation = QtWidgets.QGroupBox(Dialog_advance_setting)
        self.groupBox_receptor_preparation.setStyleSheet("font: 75 12pt \"Arial\";")
        self.groupBox_receptor_preparation.setFlat(False)
        self.groupBox_receptor_preparation.setCheckable(False)
        self.groupBox_receptor_preparation.setObjectName("groupBox_receptor_preparation")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_receptor_preparation)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.radioButton_rec_ad4 = QtWidgets.QRadioButton(self.groupBox_receptor_preparation)
        self.radioButton_rec_ad4.setChecked(True)
        self.radioButton_rec_ad4.setObjectName("radioButton_rec_ad4")
        self.gridLayout_4.addWidget(self.radioButton_rec_ad4, 0, 0, 1, 1)
        self.radioButton_rec_meeko = QtWidgets.QRadioButton(self.groupBox_receptor_preparation)
        self.radioButton_rec_meeko.setEnabled(False)
        self.radioButton_rec_meeko.setObjectName("radioButton_rec_meeko")
        self.gridLayout_4.addWidget(self.radioButton_rec_meeko, 0, 1, 1, 1)
        self.stackedWidget_rec_ad4_opt_parameters = QtWidgets.QStackedWidget(self.groupBox_receptor_preparation)
        self.stackedWidget_rec_ad4_opt_parameters.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.stackedWidget_rec_ad4_opt_parameters.setFrameShadow(QtWidgets.QFrame.Plain)
        self.stackedWidget_rec_ad4_opt_parameters.setLineWidth(1)
        self.stackedWidget_rec_ad4_opt_parameters.setObjectName("stackedWidget_rec_ad4_opt_parameters")
        self.page_rec_ad4 = QtWidgets.QWidget()
        self.page_rec_ad4.setObjectName("page_rec_ad4")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.page_rec_ad4)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.checkBox_rec_ad4_opt_parameters = QtWidgets.QCheckBox(self.page_rec_ad4)
        self.checkBox_rec_ad4_opt_parameters.setObjectName("checkBox_rec_ad4_opt_parameters")
        self.gridLayout_3.addWidget(self.checkBox_rec_ad4_opt_parameters, 0, 0, 1, 1)
        self.splitter_rec_ad4_A = QtWidgets.QSplitter(self.page_rec_ad4)
        self.splitter_rec_ad4_A.setEnabled(True)
        self.splitter_rec_ad4_A.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_rec_ad4_A.setObjectName("splitter_rec_ad4_A")
        self.checkBox_rec_ad4_A = QtWidgets.QCheckBox(self.splitter_rec_ad4_A)
        self.checkBox_rec_ad4_A.setEnabled(False)
        self.checkBox_rec_ad4_A.setObjectName("checkBox_rec_ad4_A")
        self.comboBox_rec_ad4_A_setting = QtWidgets.QComboBox(self.splitter_rec_ad4_A)
        self.comboBox_rec_ad4_A_setting.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comboBox_rec_ad4_A_setting.sizePolicy().hasHeightForWidth())
        self.comboBox_rec_ad4_A_setting.setSizePolicy(sizePolicy)
        self.comboBox_rec_ad4_A_setting.setEditable(False)
        self.comboBox_rec_ad4_A_setting.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToContents)
        self.comboBox_rec_ad4_A_setting.setFrame(True)
        self.comboBox_rec_ad4_A_setting.setObjectName("comboBox_rec_ad4_A_setting")
        self.comboBox_rec_ad4_A_setting.addItem("")
        self.comboBox_rec_ad4_A_setting.addItem("")
        self.comboBox_rec_ad4_A_setting.addItem("")
        self.comboBox_rec_ad4_A_setting.addItem("")
        self.comboBox_rec_ad4_A_setting.addItem("")
        self.gridLayout_3.addWidget(self.splitter_rec_ad4_A, 1, 0, 1, 2)
        self.splitter_rec_ad4_single = QtWidgets.QSplitter(self.page_rec_ad4)
        self.splitter_rec_ad4_single.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_rec_ad4_single.setObjectName("splitter_rec_ad4_single")
        self.checkBox_rec_ad4_C = QtWidgets.QCheckBox(self.splitter_rec_ad4_single)
        self.checkBox_rec_ad4_C.setEnabled(False)
        self.checkBox_rec_ad4_C.setObjectName("checkBox_rec_ad4_C")
        self.checkBox_rec_ad4_e = QtWidgets.QCheckBox(self.splitter_rec_ad4_single)
        self.checkBox_rec_ad4_e.setEnabled(False)
        self.checkBox_rec_ad4_e.setObjectName("checkBox_rec_ad4_e")
        self.checkBox_rec_ad4_w = QtWidgets.QCheckBox(self.splitter_rec_ad4_single)
        self.checkBox_rec_ad4_w.setEnabled(False)
        self.checkBox_rec_ad4_w.setObjectName("checkBox_rec_ad4_w")
        self.gridLayout_3.addWidget(self.splitter_rec_ad4_single, 2, 0, 1, 1)
        self.splitter_rec_ad4_p = QtWidgets.QSplitter(self.page_rec_ad4)
        self.splitter_rec_ad4_p.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_rec_ad4_p.setObjectName("splitter_rec_ad4_p")
        self.checkBox_rec_ad4_p = QtWidgets.QCheckBox(self.splitter_rec_ad4_p)
        self.checkBox_rec_ad4_p.setEnabled(False)
        self.checkBox_rec_ad4_p.setMaximumSize(QtCore.QSize(40, 16777215))
        self.checkBox_rec_ad4_p.setObjectName("checkBox_rec_ad4_p")
        self.lineEdit_rec_ad4_p = QtWidgets.QLineEdit(self.splitter_rec_ad4_p)
        self.lineEdit_rec_ad4_p.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_rec_ad4_p.sizePolicy().hasHeightForWidth())
        self.lineEdit_rec_ad4_p.setSizePolicy(sizePolicy)
        self.lineEdit_rec_ad4_p.setDragEnabled(False)
        self.lineEdit_rec_ad4_p.setObjectName("lineEdit_rec_ad4_p")
        self.gridLayout_3.addWidget(self.splitter_rec_ad4_p, 3, 0, 1, 1)
        self.splitter_rec_ad4_d = QtWidgets.QSplitter(self.page_rec_ad4)
        self.splitter_rec_ad4_d.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_rec_ad4_d.setObjectName("splitter_rec_ad4_d")
        self.checkBox_rec_ad4_d = QtWidgets.QCheckBox(self.splitter_rec_ad4_d)
        self.checkBox_rec_ad4_d.setEnabled(False)
        self.checkBox_rec_ad4_d.setMaximumSize(QtCore.QSize(40, 16777215))
        self.checkBox_rec_ad4_d.setObjectName("checkBox_rec_ad4_d")
        self.lineEdit_rec_ad4_d = QtWidgets.QLineEdit(self.splitter_rec_ad4_d)
        self.lineEdit_rec_ad4_d.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_rec_ad4_d.sizePolicy().hasHeightForWidth())
        self.lineEdit_rec_ad4_d.setSizePolicy(sizePolicy)
        self.lineEdit_rec_ad4_d.setClearButtonEnabled(False)
        self.lineEdit_rec_ad4_d.setObjectName("lineEdit_rec_ad4_d")
        self.gridLayout_3.addWidget(self.splitter_rec_ad4_d, 3, 1, 1, 1)
        self.groupBox_rec_ad4_U_setting = QtWidgets.QGroupBox(self.page_rec_ad4)
        self.groupBox_rec_ad4_U_setting.setEnabled(False)
        self.groupBox_rec_ad4_U_setting.setFlat(False)
        self.groupBox_rec_ad4_U_setting.setCheckable(True)
        self.groupBox_rec_ad4_U_setting.setChecked(False)
        self.groupBox_rec_ad4_U_setting.setObjectName("groupBox_rec_ad4_U_setting")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_rec_ad4_U_setting)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.checkBox_rec_ad4_U_nphs = QtWidgets.QCheckBox(self.groupBox_rec_ad4_U_setting)
        self.checkBox_rec_ad4_U_nphs.setChecked(True)
        self.checkBox_rec_ad4_U_nphs.setObjectName("checkBox_rec_ad4_U_nphs")
        self.gridLayout_2.addWidget(self.checkBox_rec_ad4_U_nphs, 0, 0, 1, 1)
        self.checkBox_rec_ad4_U_lps = QtWidgets.QCheckBox(self.groupBox_rec_ad4_U_setting)
        self.checkBox_rec_ad4_U_lps.setChecked(True)
        self.checkBox_rec_ad4_U_lps.setObjectName("checkBox_rec_ad4_U_lps")
        self.gridLayout_2.addWidget(self.checkBox_rec_ad4_U_lps, 0, 1, 1, 1)
        self.checkBox_rec_ad4_U_waters = QtWidgets.QCheckBox(self.groupBox_rec_ad4_U_setting)
        self.checkBox_rec_ad4_U_waters.setChecked(True)
        self.checkBox_rec_ad4_U_waters.setObjectName("checkBox_rec_ad4_U_waters")
        self.gridLayout_2.addWidget(self.checkBox_rec_ad4_U_waters, 0, 2, 1, 1)
        self.checkBox_rec_ad4_U_nonstdres = QtWidgets.QCheckBox(self.groupBox_rec_ad4_U_setting)
        self.checkBox_rec_ad4_U_nonstdres.setChecked(True)
        self.checkBox_rec_ad4_U_nonstdres.setObjectName("checkBox_rec_ad4_U_nonstdres")
        self.gridLayout_2.addWidget(self.checkBox_rec_ad4_U_nonstdres, 0, 3, 1, 1)
        self.gridLayout_3.addWidget(self.groupBox_rec_ad4_U_setting, 4, 0, 1, 2)
        self.stackedWidget_rec_ad4_opt_parameters.addWidget(self.page_rec_ad4)
        self.page_rec_meeko = QtWidgets.QWidget()
        self.page_rec_meeko.setObjectName("page_rec_meeko")
        self.checkBox_rec_meeko_opt_parameters = QtWidgets.QCheckBox(self.page_rec_meeko)
        self.checkBox_rec_meeko_opt_parameters.setGeometry(QtCore.QRect(20, 10, 191, 16))
        self.checkBox_rec_meeko_opt_parameters.setObjectName("checkBox_rec_meeko_opt_parameters")
        self.stackedWidget_rec_ad4_opt_parameters.addWidget(self.page_rec_meeko)
        self.gridLayout_4.addWidget(self.stackedWidget_rec_ad4_opt_parameters, 1, 0, 1, 2)
        self.gridLayout_8.addWidget(self.groupBox_receptor_preparation, 0, 0, 1, 1)

        self.retranslateUi(Dialog_advance_setting)
        self.stackedWidget_lig_ad4_opt_parameters.setCurrentIndex(0)
        self.stackedWidget_rec_ad4_opt_parameters.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog_advance_setting)

    def retranslateUi(self, Dialog_advance_setting):
        _translate = QtCore.QCoreApplication.translate
        Dialog_advance_setting.setWindowTitle(_translate("Dialog_advance_setting", "Dialog"))
        self.groupBox_ligand_preparation.setTitle(_translate("Dialog_advance_setting", "Ligands preparation"))
        self.radioButton_lig_ad4.setText(_translate("Dialog_advance_setting", "Autodock 4 (Default)"))
        self.radioButton_lig_meeko.setText(_translate("Dialog_advance_setting", "Meeko"))
        self.checkBox_lig_ad4_opt_parameters.setText(_translate("Dialog_advance_setting", "Optional parameters"))
        self.checkBox_lig_ad4_A.setToolTip(_translate("Dialog_advance_setting", "type(s) of repairs to make"))
        self.checkBox_lig_ad4_A.setText(_translate("Dialog_advance_setting", "-A"))
        self.comboBox_lig_ad4_A_setting.setItemText(0, _translate("Dialog_advance_setting", "None: do not make any repairs (Default)"))
        self.comboBox_lig_ad4_A_setting.setItemText(1, _translate("Dialog_advance_setting", "hydrogens: add hydrogens"))
        self.comboBox_lig_ad4_A_setting.setItemText(2, _translate("Dialog_advance_setting", "bonds_hydrogens: build bonds and add hydrogens"))
        self.comboBox_lig_ad4_A_setting.setItemText(3, _translate("Dialog_advance_setting", "bonds: build a single bond from each atom with no bonds to its closest neighbor"))
        self.checkBox_lig_ad4_R.setToolTip(_translate("Dialog_advance_setting", "index for root"))
        self.checkBox_lig_ad4_R.setText(_translate("Dialog_advance_setting", "-R"))
        self.lineEdit_lig_ad4_R.setPlaceholderText(_translate("Dialog_advance_setting", "eg: 15"))
        self.checkBox_lig_ad4_I.setToolTip(_translate("Dialog_advance_setting", "string of bonds to inactivate composed of zero-based atom indices \n"
"eg (5,13), (2,10) will inactivate atoms[5]-atoms[13] bond and atoms[2]-atoms[10] bond"))
        self.checkBox_lig_ad4_I.setText(_translate("Dialog_advance_setting", "-I"))
        self.lineEdit_lig_ad4_I.setPlaceholderText(_translate("Dialog_advance_setting", "eg (5,13), (2,10)"))
        self.checkBox_lig_ad4_p.setToolTip(_translate("Dialog_advance_setting", "preserve input charges on specific atom types, eg -p Zn -p Fe"))
        self.checkBox_lig_ad4_p.setText(_translate("Dialog_advance_setting", "-p"))
        self.lineEdit_lig_ad4_p.setPlaceholderText(_translate("Dialog_advance_setting", "Atom name (eg: Fe,Zn,C,N)"))
        self.checkBox_lig_ad4_d.setToolTip(_translate("Dialog_advance_setting", "file to contain receptor summary information"))
        self.checkBox_lig_ad4_d.setText(_translate("Dialog_advance_setting", "-d"))
        self.lineEdit_lig_ad4_d.setPlaceholderText(_translate("Dialog_advance_setting", "File name"))
        self.groupBox_lig_ad4_U_setting.setToolTip(_translate("Dialog_advance_setting", "cleanup type: nphs_lps, nphs, lps, \'\' (default is \'nphs_lps\')"))
        self.groupBox_lig_ad4_U_setting.setTitle(_translate("Dialog_advance_setting", "-U"))
        self.checkBox_lig_ad4_U_nphs.setToolTip(_translate("Dialog_advance_setting", "merge charges and remove non-polar hydrogens"))
        self.checkBox_lig_ad4_U_nphs.setText(_translate("Dialog_advance_setting", "nphs"))
        self.checkBox_lig_ad4_U_lps.setToolTip(_translate("Dialog_advance_setting", "merge charges and remove lone pairs"))
        self.checkBox_lig_ad4_U_lps.setText(_translate("Dialog_advance_setting", "lps"))
        self.groupBox_lig_ad4_B_setting.setToolTip(_translate("Dialog_advance_setting", "type(s) of bonds to allow to rotate\n"
"(default sets \'backbone\' rotatable and \'amide\' + \'guanidinium\' non-rotatable)"))
        self.groupBox_lig_ad4_B_setting.setTitle(_translate("Dialog_advance_setting", "-B"))
        self.checkBox_lig_ad4_B_backbone.setText(_translate("Dialog_advance_setting", "backbone"))
        self.checkBox_lig_ad4_B_amide.setText(_translate("Dialog_advance_setting", "amide"))
        self.checkBox_lig_ad4_B_guanidinium.setText(_translate("Dialog_advance_setting", "guanidinium"))
        self.checkBox_lig_ad4_C.setToolTip(_translate("Dialog_advance_setting", "preserve all input charges ie do not add new charges\n"
"(default is addition of gasteiger charges)"))
        self.checkBox_lig_ad4_C.setText(_translate("Dialog_advance_setting", "-C"))
        self.checkBox_lig_ad4_Z.setToolTip(_translate("Dialog_advance_setting", "inactivate all active torsions\n"
"(default is leave all rotatable active except amide and guanidinium)"))
        self.checkBox_lig_ad4_Z.setText(_translate("Dialog_advance_setting", "-Z"))
        self.checkBox_lig_ad4_g.setToolTip(_translate("Dialog_advance_setting", "attach all nonbonded fragments"))
        self.checkBox_lig_ad4_g.setText(_translate("Dialog_advance_setting", "-g"))
        self.checkBox_lig_ad4_s.setToolTip(_translate("Dialog_advance_setting", "attach all nonbonded singletons:\n"
"NB: sets attach all nonbonded fragments too\n"
"(default is not to do this)"))
        self.checkBox_lig_ad4_s.setText(_translate("Dialog_advance_setting", "-s"))
        self.checkBox_lig_ad4_w.setToolTip(_translate("Dialog_advance_setting", "assign each ligand atom a unique name: newname is original name plus its index(1-based)"))
        self.checkBox_lig_ad4_w.setText(_translate("Dialog_advance_setting", "-w"))
        self.checkBox_lig_ad4_F.setToolTip(_translate("Dialog_advance_setting", "check for and use largest non-bonded fragment (default is not to do this)"))
        self.checkBox_lig_ad4_F.setText(_translate("Dialog_advance_setting", "-F"))
        self.checkBox_meeko_lig_opt_parameters.setText(_translate("Dialog_advance_setting", "Optional parameters"))
        self.groupBox_reduce_onoff.setTitle(_translate("Dialog_advance_setting", "REDUCE 3.16 (hydrogen atoms adjustment for receptor)"))
        self.radioButton_FILP.setToolTip(_translate("Dialog_advance_setting", "add H and rotate and flip NQH groups"))
        self.radioButton_FILP.setText(_translate("Dialog_advance_setting", "FLIP"))
        self.radioButton_NOFLIP.setToolTip(_translate("Dialog_advance_setting", "add H and rotate groups with no NQH flips"))
        self.radioButton_NOFLIP.setText(_translate("Dialog_advance_setting", "NOFLIP"))
        self.radioButton_Trim.setToolTip(_translate("Dialog_advance_setting", "remove (rather than add) hydrogens"))
        self.radioButton_Trim.setText(_translate("Dialog_advance_setting", "Trim (Default)"))
        self.radioButton_BUILD.setToolTip(_translate("Dialog_advance_setting", "add H, including His sc NH, then rotate and flip groups\n"
"(except for pre-existing methionine methyl hydrogens)"))
        self.radioButton_BUILD.setText(_translate("Dialog_advance_setting", "BUILD"))
        self.pushButton_Save.setText(_translate("Dialog_advance_setting", "Save"))
        self.pushButton_Cancel.setText(_translate("Dialog_advance_setting", "Cancel"))
        self.groupBox_receptor_preparation.setTitle(_translate("Dialog_advance_setting", "Receptor preparation"))
        self.radioButton_rec_ad4.setText(_translate("Dialog_advance_setting", "Autodock 4 (Default)"))
        self.radioButton_rec_meeko.setText(_translate("Dialog_advance_setting", "Meeko"))
        self.checkBox_rec_ad4_opt_parameters.setText(_translate("Dialog_advance_setting", "Optional parameters"))
        self.checkBox_rec_ad4_A.setToolTip(_translate("Dialog_advance_setting", "type(s) of repairs to make"))
        self.checkBox_rec_ad4_A.setText(_translate("Dialog_advance_setting", "-A"))
        self.comboBox_rec_ad4_A_setting.setItemText(0, _translate("Dialog_advance_setting", "None: do not make any repairs (Default)"))
        self.comboBox_rec_ad4_A_setting.setItemText(1, _translate("Dialog_advance_setting", "hydrogens: add hydrogens"))
        self.comboBox_rec_ad4_A_setting.setItemText(2, _translate("Dialog_advance_setting", "checkhydrogens: add hydrogens only if there are none already"))
        self.comboBox_rec_ad4_A_setting.setItemText(3, _translate("Dialog_advance_setting", "bonds_hydrogens: build bonds and add hydrogens"))
        self.comboBox_rec_ad4_A_setting.setItemText(4, _translate("Dialog_advance_setting", "bonds: build a single bond from each atom with no bonds to its closest neighbor"))
        self.checkBox_rec_ad4_C.setToolTip(_translate("Dialog_advance_setting", "preserve all input charges ie do not add new charges\n"
"(default is addition of gasteiger charges)"))
        self.checkBox_rec_ad4_C.setText(_translate("Dialog_advance_setting", "-C"))
        self.checkBox_rec_ad4_e.setToolTip(_translate("Dialog_advance_setting", "delete every nonstd residue from any chain\n"
"any residue whose name is not in this list:\n"
"[\'CYS\',\'ILE\',\'SER\',\'VAL\',\'GLN\',\'LYS\',\'ASN\',\n"
"  \'PRO\',\'THR\',\'PHE\',\'ALA\',\'HIS\',\'GLY\',\'ASP\',\n"
"  \'LEU\', \'ARG\', \'TRP\', \'GLU\', \'TYR\',\'MET\',\n"
"  \'HID\', \'HSP\', \'HIE\', \'HIP\', \'CYX\', \'CSS\']\n"
"will be deleted from any chain."))
        self.checkBox_rec_ad4_e.setText(_translate("Dialog_advance_setting", "-e"))
        self.checkBox_rec_ad4_w.setToolTip(_translate("Dialog_advance_setting", "assign each receptor atom a unique name: newname is original name plus its index(1-based)"))
        self.checkBox_rec_ad4_w.setText(_translate("Dialog_advance_setting", "-w"))
        self.checkBox_rec_ad4_p.setToolTip(_translate("Dialog_advance_setting", "preserve input charges on specific atom types, eg -p Zn -p Fe"))
        self.checkBox_rec_ad4_p.setText(_translate("Dialog_advance_setting", "-p"))
        self.lineEdit_rec_ad4_p.setPlaceholderText(_translate("Dialog_advance_setting", "Atom name (eg: Fe,Zn,C,N)"))
        self.checkBox_rec_ad4_d.setToolTip(_translate("Dialog_advance_setting", "file to contain receptor summary information"))
        self.checkBox_rec_ad4_d.setText(_translate("Dialog_advance_setting", "-d"))
        self.lineEdit_rec_ad4_d.setPlaceholderText(_translate("Dialog_advance_setting", "File name"))
        self.groupBox_rec_ad4_U_setting.setTitle(_translate("Dialog_advance_setting", "-U"))
        self.checkBox_rec_ad4_U_nphs.setToolTip(_translate("Dialog_advance_setting", "merge charges and remove non-polar hydrogens"))
        self.checkBox_rec_ad4_U_nphs.setText(_translate("Dialog_advance_setting", "nphs"))
        self.checkBox_rec_ad4_U_lps.setToolTip(_translate("Dialog_advance_setting", "merge charges and remove lone pairs"))
        self.checkBox_rec_ad4_U_lps.setText(_translate("Dialog_advance_setting", "lps"))
        self.checkBox_rec_ad4_U_waters.setToolTip(_translate("Dialog_advance_setting", "remove water residues"))
        self.checkBox_rec_ad4_U_waters.setText(_translate("Dialog_advance_setting", "waters"))
        self.checkBox_rec_ad4_U_nonstdres.setToolTip(_translate("Dialog_advance_setting", "remove chains composed entirely of residues of\n"
"types other than the standard 20 amino acids"))
        self.checkBox_rec_ad4_U_nonstdres.setText(_translate("Dialog_advance_setting", "nonstdres"))
        self.checkBox_rec_meeko_opt_parameters.setText(_translate("Dialog_advance_setting", "Optional parameters"))
