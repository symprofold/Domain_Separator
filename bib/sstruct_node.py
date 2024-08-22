import ctl


def add_sstruct_to_nodes(residue_ranges, nodes):
    '''
    Integrate new secondary structure (residue_ranges) to list of existing
    nodes.
    The secondary structure is either added to its associated node, or a new
    node is created if neccessary.
    '''
    nodes_added = []
    added = False

    for n in nodes:
        n_, added_ = add_sstruct_to_node(residue_ranges, n)
        nodes_added.append(n_)

        if added_:
            added = True

    nodes = nodes_added
    check_nodes(nodes)

    if not added:
        residues_set_all = set() # set of all residues
        residues_chains = [] # residue ranges

        for residues in residue_ranges:
            residues_set_all.update(set(residues))
            residues_chains.append(tuple(residues))

        nodes.append([residues_set_all, \
                      consolidate_sstruct_ranges(residues_chains)])
            # [set of all residues, residue ranges]

    nodes = consolidate_nodes(nodes)
    nodes = sort_nodes(nodes)
    check_nodes(nodes)

    return nodes


def add_sstruct_to_node(residue_ranges0, n):
    '''
    If the secondary structure can be crosslinked with node n, add it to n.
    If it can not be crosslinked, return added = False.
    '''
    residue_ranges = consolidate_sstruct_ranges(residue_ranges0)
    added = False
    
    for residues in residue_ranges:
        residues = tuple(residues)

        intersection = set(residues) & set(n[0])

        if len(intersection) == 0:
            continue

        for residues_ in residue_ranges:
            residues_ = tuple(residues_)

            if residues != residues_:

                n[0].update(set(residues_))

                if residues_ in n[1]:
                    pass
                elif residues_ not in n[1]:
                    n[1].append(residues_)

                added = True

    if added:
        n1 = [n[0], []]
        n1[1] = consolidate_sstruct_ranges(n[1])
    else:
        n1 = [n[0], n[1]]
      
    check_sstruct_ranges(n1)
    check_sstruct_ranges1(n1[1])

    return n1, added


def consolidate_nodes(nodes):
    '''
    Merge nodes that have common residues.
    Input: list of nodes
    '''
    for i,n0 in enumerate(nodes):
        for j,n1 in enumerate(nodes):

            if j >= i:
                continue
            
            intersection = set(n0[0]) & set(n1[0])

            if len(intersection) > 0:
                n0_new, added = add_sstruct_to_node(n1[1], n0)
                nodes[i] = n0_new
                del nodes[j]

                return consolidate_nodes(nodes)

    return nodes


def consolidate_sstruct_ranges(n):
    '''
    Merge ranges that have common residues.
    Input: node
    '''
    n1 = []
    repeat = False

    for i,range0 in enumerate(n):
        toadd = True
        range0_set = set(range0)

        for j,range1 in enumerate(n1):
            range1_set = set(range1)
            intersection = range0_set & range1_set
            union = tuple(sorted(range0_set | range1_set))

            if len(intersection) == len(range0):
                toadd = False
                break

            if len(intersection) > 0:
                n1[j] = union
                repeat = True 
                toadd = False
                break

        if toadd:
            n1.append(range0)

    n1 = sort_ranges(n1)

    if repeat:
        return consolidate_sstruct_ranges(n1)

    check_sstruct_ranges1(n1)

    return n1


# Sorting
# -------

def sort_nodes(nodes):
    ''' Sort list of nodes. '''

    nodes_sort = sorted(nodes, key=lambda x: x[1][0][0])

    return nodes_sort


def sort_ranges(ranges):
    ''' Sort list of ranges. '''

    ranges_sort = sorted(ranges, key=lambda x: x[0])

    return ranges_sort


# Checks of nodes data structure
# ------------------------------

def check_nodes(nodes):
    '''
    Check data structure of nodes list.
    '''
    for n in nodes:
        if not check_sstruct_ranges(n):
            ctl.e('check_nodes: check_sstruct_ranges')
        if not check_sstruct_ranges1(n[1]):
            ctl.e('check_nodes: check_sstruct_ranges1')

    return


def check_sstruct_ranges(n):
    '''
    Check data structure of node, part: set of all residues.
    '''
    l0 = len(n[0])
    l1 = 0

    for l in n[1]:
        l1 += len(l)

    if l0 != l1:
        ctl.e(n)
        ctl.e(l0)
        ctl.e(l1)
        ctl.error('check_sstruct_ranges: check failed')

    return True


def check_sstruct_ranges1(n1):
    '''
    Check data structure of node, part: residue ranges
    '''
    set_ = set()
    l0 = 0

    for l in n1:
        l0 += len(l)
        set_ = set_ | set(l)

    l1 = len(set_)

    if l0 != l1:
        ctl.e(n1)
        ctl.e(l0)
        ctl.e(l1)
        ctl.error('check_sstruct_ranges1: check failed')

    return True


# Determination of node properties
# --------------------------------

def get_node_params(s, seq_len):
    '''
    Get parameter of node.
    '''
    order = len(s[1])
    d = get_distances(s[1])
    range_ = [s[1][0][0], s[1][-1][-1]]

    return order, range_, d


def get_node_info(s, seq_len):
    '''
    Get infotext of node.
    '''
    order, range_, d = get_node_params(s, seq_len)
    infotxt = 'order: '+str(order)+", "+str(range_[0])+'-'+str(range_[1])

    return infotxt


def nodes_on_chain(sstruct_nodes, seq_len):
    '''
    Get node positions on residue chain.
    '''
    nodes_on_chain = [-1 for i in range(0, seq_len)]
        # start from 0

    for i,n in enumerate(sstruct_nodes):
        for j,r in enumerate(n[0]):
            nodes_on_chain[r-1] = i

    return nodes_on_chain


# Determination of node properties: edges
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_distances(s):
    '''
    Get distances/sequence lengths between nodes.
    '''
    dist = []

    for i,range_ in enumerate(s):
        if i == 0:
            continue
        
        dist.append(s[i][0]-s[i-1][-1])

    return dist


def get_connections(sstruct_nodes, seq_len):
    '''
    Get connections between nodes.
    '''
    conn = []
    nodes_on_ch = nodes_on_chain(sstruct_nodes, seq_len)
    
    conn_active = False
    conn_len_count = 0

    for i,n in enumerate(nodes_on_ch):

        if n == -1:
            conn_len_count += 1

        if n != -1 and not conn_active:
            if len(conn) > 0:
                conn[-1][1] = n
                conn[-1][2] = conn_len_count

                # delete self conncetion
                if conn[-1][0] == conn[-1][1]:
                    conn.pop()

            conn.append([n, -1, -1])
            conn_len_count = 0 # reset counter
            conn_active = True
            continue
        elif n != -1 and conn_active:
            continue
        elif n == -1 and conn_active:
            conn_active = False
            continue

    if len(conn) > 0:
        conn.pop()

    return conn
