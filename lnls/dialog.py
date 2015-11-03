import warnings
warnings.filterwarnings('ignore')
import os as _os
from guidata.qt.QtGui import QApplication, QFileDialog, QWidget, QDesktopWidget
from guidata.qt.QtGui import QListView, QTreeView, QFileSystemModel, QGridLayout
from guidata.qt.QtGui import QPushButton, QLabel, QLineEdit
from guidata.qt.QtGui import QAbstractItemView
from guidata.qt.QtCore import QCoreApplication

CURDIR = _os.path.abspath(_os.path.curdir)

def directories_dialog(path=None,name='Select Directories'):
    ok = True
    def _pressed_cancel():
        nonlocal ok
        Fi.close()
        ok &= False

    path = path or CURDIR

    try:
        app = QApplication([])
    except RuntimeError:
        pass

    Fi = QFileDialog()
    Fi.setWindowTitle(name)
    Fi.setOption(Fi.DontUseNativeDialog, True)
    qr = Fi.frameGeometry()
    cp = QDesktopWidget().availableGeometry().center()
    qr.moveCenter(cp)
    Fi.move(qr.topLeft())

    Fi.setFileMode(Fi.DirectoryOnly)
    Fi.setDirectory(path)
    for view in Fi.findChildren(QListView):
        if isinstance(view.model(), QFileSystemModel):
             view.setSelectionMode(QAbstractItemView.MultiSelection)
    for view in Fi.findChildren(QTreeView):
        if isinstance(view.model(), QFileSystemModel):
             view.setSelectionMode(QAbstractItemView.MultiSelection)
    for view in Fi.findChildren(QPushButton):
        if view.text().lower().startswith('cancel'):
            view.clicked.connect(_pressed_cancel)

    Fi.show()
    QCoreApplication.instance().exec_()

    # The folder selection is also selecting its parent:
    sel_files = Fi.selectedFiles()
    sel_files2 = set(sel_files)
    for fi1 in sel_files:
        for fi2 in sel_files:
            if fi2 != fi1 and fi1 in fi2:
                sel_files2 -= {fi1}
                break

    return ok, list(sel_files2)

def input_dialog(prompt,def_answer=None,name='Type Parameters'):
    ok = False
    def _pressed_ok():
        nonlocal ok
        w.close()
        ok |= True
    def _pressed_cancel():
        nonlocal ok
        w.close()
        ok &= False

    if isinstance(prompt,str): prompt = [prompt]
    if def_answer is None: def_answer = len(prompt)*['']
    if isinstance(def_answer,str): def_answer = [def_answer]
    if len(prompt) != len(def_answer):
        raise IndexError("'prompt' and 'def_answer' must be the same length.")

    try:
        app = QApplication([])
    except RuntimeError:
        pass

    w = QWidget()
    w.setWindowTitle(name)
    grid = QGridLayout()
    grid.setSpacing(10)
    edit = []
    for i in range(len(prompt)):
        title  = QLabel(prompt[i])
        edit  += [QLineEdit()]
        if def_answer is not None: edit[i].setText(def_answer[i])
        grid.addWidget(title, 2*i,  0,1,2)# title, row,col,spanrow,spancol
        grid.addWidget(edit[i], 2*i+1, 0,1,2)
    #Ok Button
    qbtn = QPushButton('Ok', w)
    qbtn.clicked.connect(_pressed_ok)
    qbtn.resize(qbtn.sizeHint())
    grid.addWidget(qbtn, 2*(i+1), 0)
    #Cancel Button
    qbtn = QPushButton('Cancel', w)
    qbtn.clicked.connect(_pressed_cancel)
    qbtn.resize(qbtn.sizeHint())
    grid.addWidget(qbtn, 2*(i+1), 1)

    #Defining the layout of the window:
    w.setLayout(grid)
    w.resize(50, i*50)
    qr = w.frameGeometry()
    cp = QDesktopWidget().availableGeometry().center()
    qr.moveCenter(cp)
    w.move(qr.topLeft())
    w.show()
    QCoreApplication.instance().exec_()

    text = []
    for ed in edit:
        text += [ed.text()]
    return ok, text
