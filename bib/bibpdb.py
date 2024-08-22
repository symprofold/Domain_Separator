import ctl
import filesystem


'''
Module providing functions for handling pdb files.
'''

def get_f(l, field):
    '''
    Get field in line of pdb file.

    fields: 'name', 'res', 'chainid', 'resnr', 'x', 'atnr', 'element', 'bfact'
    '''
    if field == 'name':
        return l[13:16].strip()
    if field == 'res':
        return l[17:20].strip()
    if field == 'chainid':
        return l[21:22]
    if field == 'resnr':
        return int(l[22:30])
    if field == 'x':
        return float(l[30:38])
    if field == 'y':
        return float(l[38:46])
    if field == 'z':
        return float(l[46:54])
    if field == 'atnr':
        return int(l[4:11])
    if field == 'element':
        return l[77:78]
    if field == 'bfact':
        return float(l[60:66])

    return false


def set_f(l, field, nr):
    '''
    Set field in line of pdb file.

    fields: 'chainid', 'resnr', 'bfact'
    '''
    r = field

    if field == 'chainid':
        if len(nr) != 1:
            ctl.error('ERROR: chainid length')

        r = l[0:21]+str(nr)+l[22:]
    if field == 'resnr':
        r = l[0:22]+('          '+str(nr))[-4:]+l[26:]
    if field == 'bfact':
        if nr > 999.99:
            nr = 999.99
        r = l[0:60]+('          '+str(round(nr,2)))[-6:]+l[66:]

    return r


def clean_pdb(file, fileout=''):
    '''
    Clean/reformat pdb file.

    - Reduce file to the entries 'ATOM', 'TER' and 'END'.
    - Renumber chains starting with 1 and ensure continuous residue numbers.
    '''
    coord0 = filesystem.get_file(file)
    if fileout == '':
        fileout = file

    # replace 'ENDMDL' code with 'TER' to ensure the presence of a
    # termination code
    coord0 = [i.replace('ENDMDL', 'TER') for i in coord0]


    # remove duplicates of 'TER' in 2 consecutive lines
    coord = []
    prev_line = ''

    for l in coord0:
        if len(l) >= 3:
            if l[0:3] == 'TER' == (prev_line+'   ')[0:3]:
                prev_line = l
                continue

        coord.append(l)
        prev_line = l


    # remove all lines not starting with 'ATOM', 'TER' or 'END'
    coord2 = []

    for l in coord:
        if l[0:4] == 'ATOM' or l[0:3] == 'TER' or l[0:3] == 'END':
            coord2.append(l)

    coord = coord2


    # renumber chains starting with 1 and ensure continuous residue numbers
    res_i = -1 # current residue number
    res_i_renumb = 0 # current residue number after renumbering
    res_i_prev = -1 # residue number of previous line

    for i,l in enumerate(coord):
        if l[0:4] == 'ATOM':
            res_i = get_f(l, 'resnr')

            if res_i != res_i_prev:
                res_i_prev = res_i
                res_i_renumb += 1

            if res_i != res_i_renumb:
                coord[i] = set_f(l, 'resnr', res_i_renumb)

        elif l[0:3] == 'TER' or l[0:3] == 'END':
            res_i_prev = -1
            res_i_renumb = 0


    # write lines to file
    f = open(fileout, 'w')
    for l in coord:
            f.write(l) 

    f.close()

    return
