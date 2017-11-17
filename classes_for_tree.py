###########actual plotting#######################################
import info_processing
#from ivy import * 
#from ivy.interactive import *
import numpy as np
import pyqtgraph as pg
from pyqtgraph import functions as fn
from pyqtgraph import debug as debug
from pyqtgraph.Qt import QtCore, QtGui

from pyqtgraph import getConfigOption
import math
from pyqtgraph.Point import Point
from pyqtgraph.python2_3 import asUnicode

import ivy
### device coords: y incs downwards, x incs rightwards

############# custom subclasses ###########################################
class Tree_index_age_check():
    '''

    '''
    def __init__(self, start=0):
        self.index = start
        self.aged = True
        self.branched = True

    def index(self):
        return self.index
    def increment(self):
        self.index += 1
    def is_aged(self):
        return self.aged
    def set_as_non_aged(self):
        self.aged = False
    def is_branched(self):
        return self.branched
    def set_as_non_branched(self):
        self.branched = False

def process_tree_info(node, tree_info=Tree_index_age_check()):
    '''
    assign  node indexes too all nodes 
    and checks to see whether all the nodes are aged
    calculates number of descendants arising from each internal node
    '''
    node.ni = tree_info.index
    if tree_info.is_aged():
        if node.age is None:
            tree_info.set_as_non_aged()
    if tree_info.is_branched() and node.ni != 0:
        if node.length is None:
            tree_info.set_as_non_branched()

    tree_info.increment()
    if node.isleaf:
        node.meta['ndescendants'] = 0
        return 0
    else:
        desc_count = 0
        for child in node.children:
            desc_count += 1 + process_tree_info(child, tree_info)
        node.meta['ndescendants'] = desc_count
        return desc_count

def sample():
    A = ivy.tree.Node(ni=0, isroot=True, age=1000)
    B = ivy.tree.Node(ni=1, age=500)
    C = ivy.tree.Node(ni=2, age=300)
    D = ivy.tree.Node(ni=3, isleaf=True, age=0)
    E = ivy.tree.Node(ni=4, isleaf=True, age=0)
    F = ivy.tree.Node(ni=5, isleaf=True, age=250)
    G = ivy.tree.Node(ni=6, age=750)
    H = ivy.tree.Node(ni=7, isleaf=True, age=400)
    I = ivy.tree.Node(ni=8, isleaf=True, age=400)
    A.add_child(B)
    A.add_child(G)
    B.add_child(C)
    B.add_child(F)
    C.add_child(D)
    C.add_child(E)
    G.add_child(H)
    G.add_child(I)
    return A

class Window(pg.GraphicsWindow):
    def __init__(self, title, tree_file='examples/primates.newick',*args, **kwds):
        pg.GraphicsWindow.__init__(self, *args, **kwds)
        self.resize(750,1000)
        self.is_label = False
        self.tree = ivy.tree.read(tree_file)
        if not self.tree:
            self.tree = sample()
        self.tree_info = Tree_index_age_check()
        process_tree_info(self.tree, self.tree_info) 
        
        [self.positions, 
        self.labels, 
        self.node_types, 
        self.connections] = Tree_01x.setup_plot(tree=self.tree)
        self.view = VBox(window=self)

        self.axis = pg.AxisItem(orientation='left', linkView=self.view)
        self.tempaxis = pg.AxisItem(orientation='bottom', linkView=self.view)

        self.graph = Graph(window=self, viewbox=self.view)
        self.graph.setData(pos=self.positions, 
                           adj=self.connections, 
                           text=self.labels,
                           symbol='o', 
                           axisItems={'left': self.axis}) #'top': self.leafaxis
        self.view.graph = self.graph
        self.view.addItem(self.graph)
       
        self.addItem(self.axis, row=1, col=0, rowspan=6)
        self.addItem(self.view, row=1, col=1, rowspan=6, colspan=4)
        self.addItem(self.tempaxis, row=8, col=1, colspan=4)
       
        proxy_canopy = QtGui.QGraphicsProxyWidget()
        self.canopy_button = QtGui.QPushButton('Align to Canopy')
        self.canopy_button.clicked.connect(self.set_canopy_aligned)
        proxy_canopy.setWidget(self.canopy_button)
        self.addItem(proxy_canopy, row=0, col=5)

        proxy_perspective = QtGui.QGraphicsProxyWidget()
        self.pers_button = QtGui.QPushButton('Use Normal Zoom Center')
        self.pers_button.clicked.connect(self.toggle_perspective)
        proxy_perspective.setWidget(self.pers_button)
        self.addItem(proxy_perspective, row=1, col=5)
        
        proxy_labels = QtGui.QGraphicsProxyWidget()
        self.label_button = QtGui.QPushButton('Hide Node Labels')
        self.label_button.clicked.connect(self.show_labels)
        proxy_labels.setWidget(self.label_button)
        self.addItem(proxy_labels, row=2, col=5)

        branch_age_row = 3
        if self.tree_info.is_aged():
            proxy_age = QtGui.QGraphicsProxyWidget()
            self.age_button = QtGui.QPushButton('Set Node Ages')
            self.age_button.clicked.connect(self.set_aged)
            proxy_age.setWidget(self.age_button)
            self.addItem(proxy_age, row=branch_age_row, col=5) 

        if self.tree_info.is_branched():
            proxy_branch = QtGui.QGraphicsProxyWidget()
            self.branch_button = QtGui.QPushButton('Set Branch Lengths')
            self.branch_button.clicked.connect(self.set_branched)
            proxy_branch.setWidget(self.branch_button)
            if self.tree_info.is_aged():
                branch_age_row += 1               
            self.addItem(proxy_branch, row=branch_age_row, col=5)
        
        self.is_label = True
               
    def toggle_perspective(self):
        if self.view.is_canopy_pespective():
            self.view.set_canopy_pespective(False)
            self.pers_button.setText("Use Canopy as Zoom Center")
        else:
            self.view.set_canopy_pespective(True)
            self.pers_button.setText("Use Normal Zoom Center")
       
    def show_labels(self):
        current_range = self.view.viewRange()
        if self.is_label: #already showing text - want to hide
            self.label_button.setText("Show Node Labels")
            for item in self.graph.textItems:
                item.hide()
            
            if self.view.is_all_visible(self.view.viewRange()):
                self.view.autoRange()
            else:
                if self.view.is_canopy_pespective:
                    max_graph = self.graph.dataBounds(1)[1]
                    self.view.setRange(xRange=(current_range[0][0], current_range[0][1]),
                                     yRange=(current_range[1][0], max_graph))
            self.is_label = False
            
        else: #want to show labels
            self.label_button.setText("Hide Node Labels")
            for item in self.graph.textItems:
                item.show()
            
            if self.view.is_all_visible(current_range):
                self.view.autoRange()            
            else:
                if self.view.is_canopy_pespective:
                    max_view = self.view.childrenBounds()[1][1]
                    self.view.setRange(xRange=(current_range[0][0], current_range[0][1]),
                                     yRange=(current_range[1][0], max_view))
            self.is_label = True

    def update_tree_plot(self):
        '''
        could not just update the graph, had to initialize
        a new and and add it to the view to work around a bug
        that was shrinking the tree with every update

        also autorange is not correctly calculating the 
        bounds of textitem
        multiple right clicks can fix the problem
        decided to set range to max and min node coordinates instead
        user should then right click
        '''
        self.view.clear() # remove graph and text from view
        self.view.x_coords = set([coord[0] for coord in self.positions])
        self.view.y_coords = set([coord[1] for coord in self.positions])
       
        self.graph = Graph(window=self, viewbox=self.view)
        self.graph.setData(pos=self.positions, 
                           adj=self.connections, 
                           text=self.labels,                           
                           symbol='o', 
                           axisItems={'left': self.axis})
        
        self.view.addItem(self.graph)
        self.view.graph = self.graph
        
        # this avoids the textitem databounds bug
        self.view.setRange(xRange=(min(self.view.x_coords), max(self.view.x_coords)),
                           yRange=(min(self.view.y_coords), max(self.view.y_coords)))

        self.is_label = True
        self.label_button.setText("Hide Node Labels")
        
    def set_canopy_aligned(self):
        [self.positions, 
        self.labels, 
        self.node_types, 
        self.connections] = Tree_01x.setup_plot(tree=self.tree)
        self.update_tree_plot()
        
    def set_aged(self):
        [self.positions, 
        self.labels, 
        self.node_types, 
        self.connections] = Tree_01x.setup_plot(tree=self.tree, aged=True)
        self.update_tree_plot()

    def set_branched(self):
        [self.positions, 
        self.labels, 
        self.node_types, 
        self.connections] = Tree_01x.setup_plot(tree=self.tree, branched=True)
        self.update_tree_plot()


class VBox(pg.ViewBox):
    '''
    current problems:
    when aspect locked - wheel, rightdrag work /  but leftdrag jumps view and zooms
    not locked - wheel and leftdrag work 
               - right drag works when same scaling for  x and y
    '''
    def __init__(self, window, canopy_pespective=True, 
             *args, **kwds):

        pg.ViewBox.__init__(self, *args, **kwds)
        self.window = window
        self.canopy_pespective = canopy_pespective
        self.x_coords = set([coord[0] for coord in self.window.positions])
        self.y_coords = set([coord[1] for coord in self.window.positions])
        self.graph = None
        self.setBackgroundColor((255, 255, 255, 255))
        self.setMouseMode(self.PanMode)
        self.setMouseEnabled(x=True, y=True)
        self.enableAutoRange()

    def is_canopy_pespective(self):
        return self.canopy_pespective

    def set_canopy_pespective(self, boolean):
        self.canopy_pespective = boolean

    def is_all_visible(self, current_range=None):
        '''
        determine if all nodes are visible currently
        '''
        
        [[xmin, xmax], [ymin, ymax]] = current_range
        if xmin <= min(self.x_coords) and xmax >= max(self.x_coords) \
            and ymin <= min(self.y_coords) and ymax >= max(self.y_coords):
            return True
        return False
      
    def is_1_or_less(self, current_view):
        '''
        determine if one or less nodes are visible currently
        '''
        #for some reason in pyqtgraph, the viewRect's top gives lowest y visible
        # and bottom gives highest y visible 
        node_count = 0
        for [x, y] in self.window.positions:
            (xtruth, ytruth) = (False, False)

            if current_view.left() <= x and x <= current_view.right():
                xtruth = True
            if current_view.top() <= y and y <= current_view.bottom():
                ytruth = True

            if xtruth and ytruth:
                node_count += 1
                if node_count > 1:
                    return False
        return True
   
    def mouseDragEvent(self, ev, axis=None):
        ## Scale or translate based on mouse button
        ##leftdrag = pan || rightdrag = zoom

        if ev.button() & (QtCore.Qt.LeftButton | QtCore.Qt.MidButton):
            ev.accept() 
            pos = ev.pos() #comes in terms of viewport position (topleft is 0, 0)
            lastPos = ev.lastPos()
            dif = self.mapToView(pos) - self.mapToView(lastPos) #converts to graph's coord
            
            #print(self.mapFromDevice((pos-lastPos)))
            self._resetTarget()
            if self.canopy_pespective:
                self.translateBy(x=-dif.x()) 
            else:
                self.translateBy(x=-dif.x(), y=-dif.y()) 
        
        elif ev.button() & QtCore.Qt.RightButton:
            ev.accept()
            dif = ev.screenPos() - ev.lastScreenPos()
            dif = np.array([dif.x(), dif.y()])
            dif[0] *= -1
            s = 1.02 ** dif
            current_view = self.viewRect()

            #if only one node visible and wants to zoom in - not allowed
            if not(self.is_1_or_less(current_view) and s[1]<1):   

                #if all visible and wants to zoom out - not allowed
                if self.is_all_visible(current_range=self.viewRange()) and s>1:
                    pass
                else:  
                    # center is calculated in graph coords
                    tr = self.childGroup.transform()
                    tr = fn.invertQTransform(tr)
                    center = pg.QtCore.QPointF(tr.map(ev.buttonDownPos(QtCore.Qt.RightButton))) 
                        
                    if self.canopy_pespective: #changing center to ensure fixed canopy
                        # max y-value of graph aka canopy
                        max_y_coord = max(self.y_coords)
                        center = (center.x(), max_y_coord)
                    
                    self._resetTarget()
                    if (ev.modifiers() & QtCore.Qt.ShiftModifier):
                        self.scaleBy(x=s[1], center=center)
                    elif (ev.modifiers() & QtCore.Qt.ControlModifier):
                        self.scaleBy(y=s[1], center=center)
                    else:
                        #using y displacment for calculate both scales
                        self.scaleBy(x=s[1], y=s[1], center=center) 

    def wheelEvent(self, ev, axis=None):
        ev.accept()
        # s is actual scaling factor
        s = 1.02 ** (ev.delta() * self.state['wheelScaleFactor']) 
        current_view = self.viewRect()

        #if only one node visible and wants to zoom in - not allowed
        if not(self.is_1_or_less(current_view) and s<1):           
                        
            #if all visible and wants to zoom out - not allowed
            if self.is_all_visible(current_range=self.viewRange()) and s>1:
                pass
            else:
                center = pg.QtCore.QPointF(fn.invertQTransform(self.childGroup.transform()).map(ev.pos()))
                
                if self.canopy_pespective: #changing center to ensure fixed canopy
                    # max y-value of graph aka canopy
                    max_y_coord = max(self.y_coords)
                    center = pg.QtCore.QPointF(center.x(), max_y_coord)
                 
                self._resetTarget()
                if (ev.modifiers() & QtCore.Qt.ShiftModifier):
                    self.scaleBy(x=s, center=center)
                elif (ev.modifiers() & QtCore.Qt.ControlModifier):
                    self.scaleBy(y=s, center=center)
                else:
                    self.scaleBy(s=s, center=center)

    def mouseClickEvent(self, ev):
        if ev.button() == QtCore.Qt.RightButton:
            ev.accept()
            self.autoRange() # zooms out to whole tree visible
        elif ev.button() == QtCore.Qt.LeftButton:
            ev.accept()
            pos = ev.pos()
            in_pixel = self.mapToDevice(pos)
            print("pos ", pos)
          
            print("view ", self.mapToView(pos))
            print("scene ", ev.scenePos())
       
class Graph(pg.GraphItem):

    def __init__(self, window, viewbox):
        self.textItems = [] #maintain this order
        pg.GraphItem.__init__(self)
        #self.scatter.sigClicked.connect(self.clicked)
        self.view = viewbox
        self.window = window
    
    def setData(self, **kwds):
        self.text = kwds.pop('text', [])
        self.data = kwds
        if 'pos' in self.data:
            npts = self.data['pos'].shape[0]
            self.data['data'] = np.empty(npts, dtype=[('index', int)])
            self.data['data']['index'] = np.arange(npts)
        self.setTexts(self.text)
        self.updateGraph()
        
    def setTexts(self, text):
        for i in self.textItems:
            i.scene().removeItem(i)
        self.textItems = []
        for i, t in enumerate(self.text):
            if self.window.node_types[i] == 2: #if leaf-label
                angle = 60
                color = (100, 200, 100)
                anchor = (0,0.5)             
            else: #if root or internal node-label
                angle = 0
                color = (228,101,101) 
                anchor = (0,0)
            item = pg.TextItem(text=t, color=color, angle=angle, 
                          anchor=anchor)
            self.textItems.append(item)
            self.view.addItem(item)
        
    def updateGraph(self):
        pg.GraphItem.setData(self, **self.data)
        for i,item in enumerate(self.textItems):
            item.setPos(*self.data['pos'][i])

    def generatePicture(self):
        self.picture = QtGui.QPicture()
        if self.pen is None or self.pos is None or self.adjacency is None:
            return
        p = QtGui.QPainter(self.picture)
        pen = self.pen
        if pen == 'default':
            pen = getConfigOption('foreground')
        p.setPen(fn.mkPen((165,42,42), width=3))       
        
        pts = self.pos[self.adjacency]
        pts = pts.reshape((pts.shape[0]*pts.shape[1], pts.shape[2]))
        x=pts[:,0]
        y=pts[:,1]
        path = QtGui.QPainterPath()
        path.moveTo(x[0], y[0])
        for i in range(1, y.shape[0]):
            if i%2 != 0:
                c1x = x[i-1] + 0.4 * (x[i] - x[i-1])
                c1y = y[i-1] + 0.05 * (y[i] - y[i-1])
                c2x = x[i] + 0.2 * (x[i] - x[i-1])
                c2y = y[i] - 0.5 * (y[i] - y[i-1])
                c1 = pg.QtCore.QPointF(c1x, c1y)
                c2 = pg.QtCore.QPointF(c2x, c2y)
                endPoint = pg.QtCore.QPointF(x[i], y[i])
                path.cubicTo(c1, c2, endPoint)
            else:
                path.moveTo(x[i], y[i])
        p.drawPath(path)

   

    
    
