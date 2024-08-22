import ctl
import chimerax_api
import data
import sstruct_node

from bs4 import BeautifulSoup

import numpy as np
import os
import sys


def get_dssp_from_log(file_path):
    '''
    Extract secondary structure from dssp text output of ChimeraX.
    '''
    with open(file_path, 'r', encoding='utf-8') as f:
        txt = f.read()
    
    txt_object = BeautifulSoup(txt, 'html.parser')
    txt_wohtml = txt_object.get_text()

    return txt_wohtml


def get_dssp(root_path, sess, additional_paramaters=''):
    '''
    Get secondary structure using dssp in ChimeraX.
    '''
    dssp = sess.run('log clear')
    dssp = sess.run('dssp #1 report true '+additional_paramaters)
    dssp = sess.run('log save "'+root_path+'tmplog.htm"')
    dssp_txt = get_dssp_from_log(root_path+'tmplog.htm')
    os.unlink(root_path+'tmplog.htm')

    dssp = dssp_txt.split("\n")

    return dssp


def dssp_matrix(dssp, seq_len, maxdist_betasheet, allowed_resids=False):
    '''
    Create matrices with secondary structures.

    Returns:
        dssp_matr: matrix with crosslinks between all residues of a complete
                - beta sheet (several strands) or
                - alpha helix

        dssp_matr_betasheet_interstrand: matrix with crosslinks between all
                inter strand residues of a 
                - beta sheet (several strands)
                but without crosslinks within intra strand residues

        helices: list of alpha helix ranges 

        sstruct_nodes: list of beta sheet nodes
    '''
    dssp_matr = np.zeros((seq_len, seq_len))
        # first res is 0

    dssp_matr_betasheet_interstrand = np.zeros((seq_len, seq_len))
        # first res is 0

    # build list of secondary structure nodes (beta sheets)
    sstruct_nodes = []
        # first res is 0

    helices = []
        # first res is 0

    strands_partners = []
        # list of strand partners
        # first res is 1


    for d in dssp:

        # antiparallel beta sheets
        # ~~~~~~~~~~~~~~~~~~~~~~~~
        if ' antiparallel ' in d:
            ranges = d.split(' antiparallel ')
            range0 = ranges[0].split(' -> ')
            begin0 = int(range0[0].split(':')[1].strip())
            end0 = int(range0[1].split(':')[1].strip())

            range1 = ranges[1].split(' -> ')
            begin1 = int(range1[0].split(':')[1].strip())
            end1 = int(range1[1].split(':')[1].strip())

            if allowed_resids != False:
                range0_filtered = data.filter_range( \
                    [begin0, end0], allowed_resids[0] ^ allowed_resids[1], 2)
                if range0_filtered == []:
                    continue

                begin0 = range0_filtered[0][0]
                end0 = range0_filtered[-1][1]

                range1_filtered = data.filter_range( \
                    [begin1, end1], allowed_resids[0] ^ allowed_resids[1], 2)
                if range1_filtered == []:
                    continue

                begin1 = range1_filtered[0][0]
                end1 = range1_filtered[-1][1]


            r0 = [r for r in range(begin0, end0+1)]
            r1 = [r for r in range(begin1, end1+1)]

            if allowed_resids != False:
                if (len(r0) == 0 or len(r1) == 0):
                    continue

                if (len(r0) < 3 or len(r1) < 3):
                    continue

            sstruct_nodes = sstruct_node.add_sstruct_to_nodes(\
                                [tuple(r0), tuple(r1)], sstruct_nodes)

            strands_partners.append([tuple(r0), tuple(r1)])

            for r0 in range(begin0, end0+1):
                for r0_ in range(begin0, end0+1):
                    dssp_matr[r0-1][r0_-1] = 1

                for r1 in range(begin1, end1+1):
                    dssp_matr[r0-1][r1-1] = 1
                    dssp_matr[r1-1][r0-1] = 1
                    dssp_matr_betasheet_interstrand[r0-1][r1-1] = 1
                    dssp_matr_betasheet_interstrand[r1-1][r0-1] = 1

            for r1 in range(begin1, end1+1):
                for r1_ in range(begin1, end1+1):
                    dssp_matr[r1-1][r1_-1] = 1

            # Include external chain sections between strands of a beta sheet
            # up to a maximum length of maxdist_betasheet.
            # Thus, further secondary structures can also be included.
            dist_betasheet = max(begin0, begin1)-min(end0, end1)

            if dist_betasheet <= maxdist_betasheet:

                for r0 in range(begin0, end1+1):
                    for r0_ in range(begin0, end1+1):
                        dssp_matr[r0-1][r0_-1] = 1


        # parallel beta sheets
        # ~~~~~~~~~~~~~~~~~~~~
        elif ' parallel ' in d:
            ranges = d.split(' parallel ')
            range0 = ranges[0].split(' -> ')
            begin0 = int(range0[0].split(':')[1].strip())
            end0 = int(range0[1].split(':')[1].strip())

            range1 = ranges[1].split(' -> ')
            begin1 = int(range1[0].split(':')[1].strip())
            end1 = int(range1[1].split(':')[1].strip())

            if allowed_resids != False:
                range0_filtered = data.filter_range( \
                    [begin0, end0], allowed_resids[0] ^ allowed_resids[1], 2)
                if range0_filtered == []:
                    continue

                begin0 = range0_filtered[0][0]
                end0 = range0_filtered[-1][1]

                range1_filtered = data.filter_range( \
                    [begin1, end1], allowed_resids[0] ^ allowed_resids[1], 2)
                if range1_filtered == []:
                    continue

                begin1 = range1_filtered[0][0]
                end1 = range1_filtered[-1][1]


            r0 = [r for r in range(begin0, end0+1)]
            r1 = [r for r in range(begin1, end1+1)]

            if allowed_resids != False:
                if (len(r0) == 0 or len(r1) == 0):
                    continue

                if (len(r0) < 3 or len(r1) < 3):
                    continue

            sstruct_nodes = sstruct_node.add_sstruct_to_nodes(\
                                [tuple(r0), tuple(r1)], sstruct_nodes)

            strands_partners.append([tuple(r0), tuple(r1)])

            for r0 in range(begin0, end0+1):
                for r0_ in range(begin0, end0+1):
                    dssp_matr[r0-1][r0_-1] = 1

                for r1 in range(begin1, end1+1):
                    dssp_matr[r0-1][r1-1] = 1
                    dssp_matr[r1-1][r0-1] = 1
                    dssp_matr_betasheet_interstrand[r0-1][r1-1] = 1
                    dssp_matr_betasheet_interstrand[r1-1][r0-1] = 1

            for r1 in range(begin1, end1+1):
                for r1_ in range(begin1, end1+1):
                    dssp_matr[r1-1][r1_-1] = 1

            # Include external chain sections between strands of a beta sheet
            # up to a maximum length of maxdist_betasheet.
            # Thus, further secondary structures can also be included.
            dist_betasheet = max(begin0, begin1)-min(end0, end1)

            if dist_betasheet <= maxdist_betasheet:
                for r0 in range(begin0, end1+1):
                    for r0_ in range(begin0, end1+1):
                        dssp_matr[r0-1][r0_-1] = 1


        # alpha helices
        # ~~~~~~~~~~~~~
        elif ' -> ' in d:
            rang = d.split(' -> ')
            begin = int(rang[0].split(':')[1].strip())
            end = int(rang[1].split(':')[1].strip())

            if allowed_resids != False:
                ranges_filtered = data.filter_range([begin, end], \
                                                    allowed_resids[1])
                if ranges_filtered == []:
                    continue
            else:
                ranges_filtered = [[begin, end]]

            for range_filtered in ranges_filtered:
                begin = range_filtered[0]
                end = range_filtered[1]

                for r0 in range(begin, end+1):
                    for r0_ in range(begin, end+1):
                        dssp_matr[r0-1][r0_-1] = 1

                helices.append([begin, end])

    return dssp_matr, dssp_matr_betasheet_interstrand, helices, sstruct_nodes
