import ctl
import chimerax_api
import data
import sstruct
import sstruct_node
import visual

import copy
import numpy as np


def run_sstruct_scaffold(dssp, dssp_le, seq_len, params, options, sess):
    '''
    Build secondary structure scaffold bycrosslinking of all residues inside
    beta sheets and alpha helices.
    '''
    dssp_matr, dssp_matr_betasheet_interstrand, helix_ranges, \
        sstruct_nodes = \
            sstruct.dssp_matrix(dssp, seq_len, params['maxdist_betasheet'])
        # dssp_matr, dssp_matr_betasheet_interstrand:
        # first res is 0        
        #
        # dssp_matr: matrix with crosslinks between all residues of a complete
        #         - beta sheet (several strands) or
        #         - alpha helix
        #
        # dssp_matr_betasheet_interstrand: matrix with crosslinks between all
        #         inter strand residues of a 
        #         - beta sheet (several strands)
        #         but without crosslinks within intra strand residues
        #
        # helix_ranges: list of alpha helix ranges 
        #
        # sstruct_nodes: list of beta sheet nodes
        #
        # dssp_matr is used for crosslinking of secondary structure

    dssp_matr0 = copy.deepcopy(dssp_matr)

    if options['verbous']:
        ctl.p(str(len(sstruct_nodes))+' secondary structure nodes')

        for s in sstruct_nodes:
            info = sstruct_node.get_node_info(s, seq_len)
            ctl.p(info)

    # crosslink entries in dssp_matr to secondary structure scaffold
    dssp_matr2 = crosslink_dssp_matr(dssp_matr)
    dssp_matr_betasheet_interstrand2 = crosslink_dssp_matr( \
                                            dssp_matr_betasheet_interstrand)
        # first res is 0

    # domain ranges after crosslinking
    sstruct_sect_ranges2 = data.matr_to_ranges(dssp_matr2)
        # first res is 1

    sess.run('color #1.1 red')
    visual.color_dom_sections(sstruct_sect_ranges2, 'gray', \
                              options['dom_coloring'], sess)
    if options['show_steps']:
        sess.run('combine #1 close false modelId #10')
        sess.run('rename #10 crosslink')

    sstruct_sect_ranges = data.matr_to_ranges(dssp_matr)
    sstruct_sect_ranges = sstruct_node.sort_ranges(sstruct_sect_ranges)
    sstruct_sect_set = set(data.ranges_to_resid_list(sstruct_sect_ranges))
    helix_resids = set(data.ranges_to_resid_list(helix_ranges))

    dssp_matr_le, dssp_matr_betasheet_interstrand_le, helix_ranges_le, \
        sstruct_nodes_le = \
            sstruct.dssp_matrix(dssp_le, seq_len, \
                params['maxdist_betasheet'], [sstruct_sect_set, helix_resids])

    sstruct_sect_ranges_le = data.matr_to_ranges(dssp_matr_le)
    sstruct_sect_ranges_le = sstruct_node.sort_ranges(sstruct_sect_ranges_le)

    # crosslink entries in dssp_matr_le to secondary structure scaffold
    dssp_matr2_le = crosslink_dssp_matr(dssp_matr_le)
        # first res is 0

    # domain ranges after crosslinking
    sstruct_sect_le = data.matr_to_ranges(dssp_matr2_le)
        # first res is 1

    sstruct_sect_le_set = set(data.ranges_to_resid_list(sstruct_sect_le))
    sstruct_sect_ranges = data.resids_to_ranges( \
            sstruct_sect_set | sstruct_sect_le_set)

    return dssp_matr0, dssp_matr2, sstruct_nodes, sstruct_sect_ranges, \
           helix_ranges


def crosslink_dssp_matr(dssp_matr):
    ''' Crosslink entries in dssp_matr to secondary structure scaffold. '''

    dssp_matr_sums = np.sum(dssp_matr, axis=0)
    seq_len = len(dssp_matr_sums)
    sstruct_scaffold = np.zeros((seq_len, seq_len))

    for x in range(0, seq_len):
        for y in range(0, seq_len):
            if dssp_matr_sums[x] != 0 and dssp_matr_sums[y] != 0:
                sstruct_scaffold[x][y] = 1

    return sstruct_scaffold
