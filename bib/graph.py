import networkx as nx


def separability(g):
    '''
    Determine if removal of edges lead to disconnection of the network.
    '''
    separability = {}
    separability_global = False

    for edge in g.edges():
        g_modified = g.copy()
        g_modified.remove_edge(*edge)

        if nx.is_weakly_connected(g_modified):
            separability[edge] = False
        else:
            separability[edge] = True
            separability_global = True

    return separability, separability_global
