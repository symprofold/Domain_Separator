import ctl


def get_num_betasheet(dom_range, sstruct_nodes, max_strands=1000):
    '''
    Get number of beta sheets in a given residue range.
    '''
    betasheet_num = 0
    betasheet_ids = []

    for r in range(dom_range[0], dom_range[1]+1):
        for i,betasheet in enumerate(sstruct_nodes):
            if r in betasheet[0]:
                if i not in betasheet_ids:
                    if len(betasheet[1]) <= max_strands:
                        betasheet_ids.append(i)

    betasheet_num = len(betasheet_ids)

    return betasheet_num


def get_num_betasheet_strands(dom_range, sstruct_nodes):
    '''
    Get number of beta sheet strands in a given residue range.
    '''
    strands_num = 0
    strand_ids = []

    for r in range(dom_range[0], dom_range[1]+1):
        for i,betasheet in enumerate(sstruct_nodes):
            if r in betasheet[0]:
                for j,strand in enumerate(betasheet[1]):
                    if r in strand:
                        if (i,j) not in strand_ids:
                            strand_ids.append((i,j))

    strand_ids = len(strand_ids)

    return strand_ids


def get_num_helices(dom_range, helices):
    '''
    Get number of alpha helices in a given residue range.
    '''
    helices_num = 0
    helices_ids = []

    for r in range(dom_range[0], dom_range[1]+1):
        for i,helix in enumerate(helices):
            if r in range(helix[0], helix[1]+1):
                if i not in helices_ids:
                        helices_ids.append(i)

    helices_num = len(helices_ids)
    
    return helices_num
