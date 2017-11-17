import numpy as np


#required module 'bio' - installed on anaconda prompt with
#                        command 'conda install biopython'

# r.leaves() # gives all leaf nodes Node(979, leaf, 'Poitea') in list
# r.children # gives current node r's children in list
# f = treefig(r) 
# f.ladderize() # commands to view tree in pdf

##### current changes to paper:
## replaced mult by h_N with w_N in deepest distance formula
## replaced w''_c with w''_P in h''_c equation
## replaced s_C with s_N

##### calculating relative bounding box ####################

def assign_relative_box(node, relative_box_dict={}):
    '''
    the node width is determined by average the width's
    suggested by the paper's formula where all children of a 
    parent node get equal space (called equal_share here) and 
    my own formula where the a node gets space proportional to
    the ratio of the number of a node's own descendants 
    (node.meta['ndescendants']) to the number of descendants of itself and 
    all its siblings ( called total_sibling_descendants ). This
    latter space is called desecendant_based here
    The node_height is determined by own formula where node's with more
    children get more node height
    '''
    if not node.isroot:
        equal_share = 1 / node.parent.nchildren
        total_sibling_descendants = node.parent.meta['ndescendants'] - node.parent.nchildren
        if node.isleaf:
            desecendant_based = 0
        else:
            desecendant_based = node.meta['ndescendants'] / total_sibling_descendants
        if total_sibling_descendants == 0: ## aka if all siblings are leaf nodes
            node_width = equal_share
        else:
            node_width = (equal_share + desecendant_based) / 2

    if node.isleaf:
        relative_box_dict[node.ni] = {'width': node_width,  
                                    'height': 0,
                                    'y_position': 0}
        return relative_box_dict

    if node.isroot:
        pass
    else: # an internal non-root node
        node_height = node_width * (1 + node.nchildren/5)
        relative_box_dict[node.ni] = {'width': node_width, 
                                    'height': node_height,
                                    'y_position': node_height}
        
    x_pos = 0
    for i, child in enumerate(node.children):
        assign_relative_box(child, relative_box_dict)
        relative_box_dict[child.ni]['x_position'] = x_pos
        x_pos += relative_box_dict[child.ni]['width']
        
    return relative_box_dict

 ##### top alightment to canopy ####################

def dist_to_deepest(node, relative_box_dict, distances_dict={}):
    '''
    calculates the distance between a node and the highest tip of the
    phylogeny
    '''
    if node.isleaf:
        distances_dict[node.ni] = 0
        return distances_dict
        
    D_N_set = set() # set of all distances d_c where c is child of N
    for child in node.children:
        dist_to_deepest(child, relative_box_dict, distances_dict)
        D_N_set.add(distances_dict[child.ni])
    if not node.isroot:
        node_height = relative_box_dict[node.ni]['height'] 
        node_width = relative_box_dict[node.ni]['width'] 
        #below - is it multiply by h_n or w_n? paper mistake?
        d_N = node_height + max(D_N_set) * node_width
        distances_dict[node.ni] = d_N
        return distances_dict
    else:
        distances_dict[node.ni] = float("NAN")
        return distances_dict


def shift_values(node, relative_box_dict, distances_dict, shift_dict={}):
    '''
    calculates the amount that a node's position needs to be shifted
    so that the whole tree can be aligned to the canopy
    '''
    if not node.isleaf:
        for child in node.children:
            shift_values(child, relative_box_dict, distances_dict, shift_dict)
        
    if not node.isroot:
        D_P_N_set = set()
        for sibling in node.parent.children:
            D_P_N_set.add(distances_dict[sibling.ni])
        shift = max(D_P_N_set) - distances_dict[node.ni] 
        shift_dict[node.ni] = shift
    else: #no shift needed for root
        shift_dict[node.ni] = 0

    return shift_dict

def shift_to_canopy(relative_box_dict, shift_dict):
    '''
    changes the yposition of all node's such that the whole tree
    is aligned to the canopy
    '''
    ## here is it sc or sn? could paper be wrong?
    for node in relative_box_dict:
        relative_box_dict[node]['y_position'] += shift_dict[node]
    return relative_box_dict


 ##### transformation to absolute values #######################################

def branch_length_coordinates(node, relative_box_dict, root_width,
                      absolute_dict={}, connection_data=[]):
    '''
    recalculates the y position of the nodes based on
    branch length information and produces the absolute coordinates
    for the whole tree
    node_type: 0 - root node, 1 - internal node, 2 - leaf node
    '''
    if node.isroot:
        absolute_dict[node.ni] = {'abs_width': root_width,  
                                'abs_x_box': 0,
                                'abs_x_node': root_width * 0.5,
                                'abs_y_node': 0,
                                'label': node.label,
                                'node_type': 0}

    else:   
        abs_width = absolute_dict[node.parent.ni]['abs_width'] * relative_box_dict[node.ni]['width']
        abs_x_box = absolute_dict[node.parent.ni]['abs_x_box'] + relative_box_dict[node.ni]['x_position'] * absolute_dict[node.parent.ni]['abs_width']
        abs_y_node = absolute_dict[node.parent.ni]['abs_y_node'] + node.length
        absolute_dict[node.ni] = {'abs_width': abs_width,  
                        'abs_x_box': abs_x_box,
                        'abs_x_node': abs_x_box + 0.5 * abs_width,
                        'abs_y_node': abs_y_node,
                        'label': node.label}
        if node.isleaf:
            absolute_dict[node.ni]['node_type'] = 2
            return [absolute_dict, connection_data]

        else: # internal node 
            absolute_dict[node.ni]['node_type'] = 1

    #following will be run for root and internal nodes only
    for child in node.children:
        connection_data.append([node.ni, child.ni])
        branch_length_coordinates(child, relative_box_dict, root_width,
            absolute_dict, connection_data)

    return [absolute_dict, connection_data]

def node_aged_coordinates(node, relative_box_dict, root_width,
                      absolute_dict={}, connection_data=[]):
    '''
    recalculates the y position of the nodes based on
    age of each node and produces the absolute coordinates
    for the whole tree
    node_type: 0 - root node, 1 - internal node, 2 - leaf node
    '''
    if node.isroot:
        absolute_dict[node.ni] = {'abs_width': root_width,  
                                'abs_x_box': 0,
                                'abs_x_node': root_width * 0.5,
                                'abs_y_node': -node.age,
                                'label': node.label,
                                'node_type': 0}

    else:   
        abs_width = absolute_dict[node.parent.ni]['abs_width'] * relative_box_dict[node.ni]['width']
        abs_x_box = absolute_dict[node.parent.ni]['abs_x_box'] + relative_box_dict[node.ni]['x_position'] * absolute_dict[node.parent.ni]['abs_width']

        absolute_dict[node.ni] = {'abs_width': abs_width,  
                        'abs_x_box': abs_x_box,
                        'abs_x_node': abs_x_box + 0.5 * abs_width,
                        'abs_y_node': -node.age,
                        'label': node.label}

        if node.isleaf:
            absolute_dict[node.ni]['node_type'] = 2
            return [absolute_dict, connection_data]

        else: # internal node 
            absolute_dict[node.ni]['node_type'] = 1

    #following will be run for root and internal nodes only
    for child in node.children:
        connection_data.append([node.ni, child.ni])
        node_aged_coordinates(child, relative_box_dict, root_width,
            absolute_dict, connection_data)

    return [absolute_dict, connection_data]

def non_aged_coordinates(node, relative_box_dict, root_width, root_height,
                      x_offset=0, y_offset=0, absolute_dict={}, connection_data=[]):
    '''
    produces the absolute coordinates
    for the whole tree where the tree is aligned to the canopy
    i.e. it is not aged or have branch lengths
    node_type: 0 - root node, 1 - internal node, 2 - leaf node
    '''
    if node.isroot:
        absolute_dict[node.ni] = {'abs_width': root_width, 
                                'abs_height': root_height,
                                'abs_x_box': x_offset,
                                'abs_y_box': y_offset,
                                'abs_x_node': x_offset + root_width * 0.5,
                                'abs_y_node': y_offset - root_height,
                                'label': node.label,
                                'node_type': 0}

     
    else:   
        abs_width = absolute_dict[node.parent.ni]['abs_width'] * relative_box_dict[node.ni]['width']
        # for height is it h''_c = w''_c * h'_c as paper says
        # or h''_c = w''_p * h'_c ??
        abs_height = absolute_dict[node.parent.ni]['abs_width'] * relative_box_dict[node.ni]['height']
        abs_x_box = absolute_dict[node.parent.ni]['abs_x_box'] + relative_box_dict[node.ni]['x_position'] * absolute_dict[node.parent.ni]['abs_width']
        abs_y_box = absolute_dict[node.parent.ni]['abs_y_box'] + relative_box_dict[node.ni]['y_position'] * absolute_dict[node.parent.ni]['abs_width'] 
        abs_x_node = abs_x_box + 0.5 * abs_width
        if node.isleaf:
            abs_y_node = abs_y_box
        else: 
            abs_y_node = abs_y_box - abs_height
        
        absolute_dict[node.ni] = {'abs_width': abs_width, 
                                'abs_height': abs_height,
                                'abs_x_box': abs_x_box,
                                'abs_y_box': abs_y_box,
                                'abs_x_node': abs_x_node,
                                'abs_y_node': abs_y_node,
                                'label': node.label}

        if node.isleaf:          
            absolute_dict[node.ni]['node_type'] = 2
            return [absolute_dict, connection_data]

        else: # internal node 
            absolute_dict[node.ni]['node_type'] = 1

    #following will be run for root and internal nodes only
    for child in node.children:
        connection_data.append([node.ni, child.ni])
        non_aged_coordinates(child, relative_box_dict, root_width,
            root_height, x_offset, y_offset, absolute_dict, connection_data)

    return [absolute_dict, connection_data]


### get absolute coordinates and connections ################################

def setup_plot(tree, root_width=100, root_height=100, x_offset=0, 
              y_offset=0, aged=False, branched=False):

    relative_coord = assign_relative_box(tree, {})
    if aged: #node's y-position  informed by its age
        abs_conn = node_aged_coordinates(tree, relative_coord, root_width, {}, [])
    
    elif branched: #y-position  informed by branch lengths
        abs_conn = branch_length_coordinates(tree, relative_coord, root_width, {}, [])

    else: # y-position has no true information, tree will be aligned to canopy
        distances = dist_to_deepest(tree, relative_coord, {})
        shifts = shift_values(tree, relative_coord, distances, {})
        relative_coord = shift_to_canopy(relative_coord, shifts)
        abs_conn = non_aged_coordinates(tree, relative_coord, root_width, 
                                    root_height, x_offset, y_offset, {}, [])


    positions = np.array([[val['abs_x_node'], val['abs_y_node']] for (key,val) in abs_conn[0].items()],
                dtype=float)    
    labels = [val['label'] for (key,val) in abs_conn[0].items()]
    node_types = [val['node_type'] for (key,val) in abs_conn[0].items()]
    connections = np.array(abs_conn[1])
    
    return [positions, labels, node_types, connections]



