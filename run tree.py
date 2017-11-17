
import pyqtgraph as pg
#from pyqtgraph import functions as fn
#from pyqtgraph import Point
from pyqtgraph.Qt import QtCore, QtGui
#import pyqtgraph.examples
import sys

import info_processing
import classes_for_tree 

if len(sys.argv) > 1:
	tree_file = sys.argv[1]

	title = tree_file
else:
	tree_file = None
	title = 'Tree'
	
app = QtGui.QApplication([])
pg.setConfigOptions(antialias=True) # for prettier graphs
win = classes_for_tree.Window(title=title, tree_file=tree_file)

## Start Qt event loop unless running in interactive mode or using pyside.
'''
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()

'''
app.exec_()
