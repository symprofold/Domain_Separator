import ctl
import check_domain
import chimerax_api
import data
import domain_property
import fasta
import filesystem
import graph
import linker
import merging_score
import networkx as nx
import specialdomains
import sstruct
import sstruct_node
import sstruct_scaffold
import visual

from modelreg import ModelReg

import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


def separate_domains(f, params, options, sess):
    '''
    Separate domains of the model in input file f.
    '''
    root_path = '/'.join(f.split('/')[:-1])+'/'

    # cache for already determined surface areas and linker midpoints
    cache = {'surface_area_dict': {}, 'linker_midpoints': {}}

    sess.run('close session')
    sess.run('set bgColor white')

    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    sess.run('open "'+f+'"')


    # determine secondary structure data using dssp in ChimeraX
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dssp_le = sstruct.get_dssp(root_path, sess, ' energyCutoff -0.2')
            # low energy cutoff
    dssp = sstruct.get_dssp(root_path, sess)
            # standard energy cutoff

    # set correct model id for relevant chain to (1,1) regardless of whether
    # monomer or multimer (e.g. dimer) is provided by input file
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    seq0 = sess.get_seq((1))
    seq_len0 = len(seq0)
    sess.run('split #1')

    seq1 = sess.get_seq((1,1))
    seq_len1 = len(seq1)

    if seq_len0 > 0 and seq_len1 == 0:
        sess.run('rename #1 id #1.1')

    seq = sess.get_seq((1,1))
    seq_len = len(seq) # sequence length


    # step 1:
    # crosslink all residues inside beta sheets and alpha helices
    # ===========================================================
    dssp_matr0, dssp_matr2, sstruct_nodes, sstruct_sect_ranges, \
            helix_ranges = \
        sstruct_scaffold.run_sstruct_scaffold(
            dssp, dssp_le, seq_len, params, options, sess)

    model_id = 10

    dom_ranges = copy.deepcopy(sstruct_sect_ranges)
    dom_ranges_len = len(dom_ranges)


    # step 2:
    # determine separability of nodes from beta sheet graph
    # =====================================================

    # build graph of beta sheet crosslinks
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    g = nx.MultiDiGraph()
    node_sizes = {}

    for n in range(1, len(sstruct_nodes)+1):
        order, range_, d = sstruct_node.get_node_params(sstruct_nodes[n-1], \
                                                        seq_len)
        g.add_node(n)
        node_sizes[n] = order*100

    conn = sstruct_node.get_connections(sstruct_nodes, seq_len)

    for i,c in enumerate(conn):
        g.add_edge(c[0]+1, c[1]+1)

    # plot graph
    if options['show_plot']:
        layout = nx.spring_layout(g)
        nx.draw(g, layout, node_size=[node_sizes[n] for n in g.nodes()], \
                connectionstyle='arc3, rad = 0.1', with_labels=True)
        plt.show()


    # determine separability of nodes from graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    separability, separability_global = graph.separability(g)


    # step 3:
    # merging of domain sections where both merging candidates reach the
    # minimum score
    # ==================================================================
    dom_ranges = copy.deepcopy(sstruct_sect_ranges)
    dom_ranges_len = len(dom_ranges)


    # create list merging_score_bothabove of merging candidates where both
    # merging partners reach the minimum score
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    merging_scores, merging_score_relations, merging_score_bothabove, \
        cache = \
            merging_score.get_merging_scores( \
                dom_ranges, helix_ranges, cache, params, sess)


    # merge all merging candidates in merging_score_bothabove
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    merged_resids = set(data.ranges_to_resid_list(dom_ranges))

    for rangeid,score in enumerate(merging_score_bothabove):
        if score != -1:
            merged_resids |= set(data.ranges_to_resid_list( \
                    data.merge_ranges(dom_ranges, rangeid, rangeid+1)))

    merged_ranges = data.resids_to_ranges(merged_resids)
    dom_ranges = merged_ranges


    sess.run('color #1.1 red')
    visual.color_dom_sections(dom_ranges, 'gray', \
                              options['dom_coloring'], sess)
    if options['show_steps']:
        model_id += 1
        sess.run('combine #1 close false modelId #'+str(model_id))
        sess.run('rename #'+str(model_id)+' "both minscore"')


    dom_ranges_prev_len = -1
    i_repeat = -1

    while len(dom_ranges) != dom_ranges_prev_len:
        i_repeat += 1
        dom_ranges_prev_len = len(dom_ranges)


        # step 4:
        # iterative merging of domain sections where at least one of both merging
        # candidates reaches the minimum score
        # =======================================================================
        remaining_iterations = 1000
        dom_ranges_len = len(dom_ranges)

        for iteration in range(0, len(dom_ranges)):
            ctl.p('iteration '+str(iteration))

            if remaining_iterations < 1:
                ctl.p('finished after '+str(iteration)+' iterations')
                break


            # create list merging_scores of merging candidates where at least one
            # of both merging candidates reaches the minimum score
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            merging_scores, merging_score_relations, merging_score_bothabove, \
                cache = \
                    merging_score.get_merging_scores( \
                        dom_ranges, helix_ranges, cache, params, sess)


            # merge merging candidate with the highest score in merging_scores
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            range_id_max, score_max = \
                    merging_score.get_highest_score(merging_scores)

            if score_max < 0:
                break

            merged_ranges = data.merge_ranges(dom_ranges, range_id_max, \
                                              range_id_max+1)


            # if a domain section was merged in this iteration
            if len(merged_ranges) != dom_ranges_len:
                dom_ranges_len = len(merged_ranges)
                dom_ranges = merged_ranges
                remaining_iterations = dom_ranges_len
            else:
                remaining_iterations -= 1


        sess.run('color #1.1 red')
        visual.color_dom_sections(dom_ranges, 'gray', \
                                  options['dom_coloring'], sess)
        if options['show_steps']:
            model_id += 1
            sess.run('combine #1 close false modelId #'+str(model_id))
            sess.run('rename #'+str(model_id)+' iteration')


        # step 5: triple merging
        # iterative merging of three adjacent domain sections:
        # Two domain sections, separated by a third domain section in the
        # middle are merged, when at least one of the two domain sections/
        # merging candidates reaches the minimum score
        # ================================================================
        remaining_iterations = 1000
        dom_ranges_len = len(dom_ranges)

        for iteration in range(0, dom_ranges_len**2):
            ctl.p('iteration '+str(iteration))

            if remaining_iterations < 1:
                ctl.p('finished after '+str(iteration)+' iterations')
                break


            # create list merging_scores of merging candidates where at least
            # one of both merging candidates reaches the minimum score
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            merging_scores, merging_score_relations, \
                merging_score_bothabove, cache = \
                    merging_score.get_merging_scores( \
                        dom_ranges, helix_ranges, cache, params, sess, 1)


            # merge merging candidate with the highest score in merging_scores
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            range_id_max, score_max = \
                    merging_score.get_highest_score(merging_scores)

            if score_max < 0:
                break

            merged_ranges = data.merge_ranges(dom_ranges, range_id_max, \
                                              range_id_max+1+1, 1)


            # if a domain section was merged in this iteration
            if len(merged_ranges) != dom_ranges_len:
                dom_ranges_len = len(merged_ranges)
                dom_ranges = merged_ranges
                remaining_iterations = dom_ranges_len
            else:
                remaining_iterations -= 1


        sess.run('color #1.1 red')
        visual.color_dom_sections(dom_ranges, 'gray', \
                                  options['dom_coloring'], sess)
        if options['show_steps']:
            model_id += 1
            sess.run('combine #1 close false modelId #'+str(model_id))
            sess.run('rename #'+str(model_id)+' "iteration gap"')


        # step 6: triple merging
        # iterative merging of three adjacent domain sections:
        # Two adjacent domain sections are merged with a third, when either
        # the two domain sections taken together as one merging candidate or
        # the third domain as the other merging candidate reaches the
        # minimum score
        # ==================================================================
        remaining_iterations = 1000
        dom_ranges_len = len(dom_ranges)

        for iteration in range(0, len(dom_ranges)):
            ctl.p('iteration '+str(iteration))

            if remaining_iterations < 1:
                ctl.p('finished after '+str(iteration)+' iterations')
                break

            for merge_range_id in range(0, len(dom_ranges)-2):
                strands0 = domain_property.get_num_betasheet_strands( \
                        dom_ranges[merge_range_id], sstruct_nodes)
                strands1 = domain_property.get_num_betasheet_strands( \
                        dom_ranges[merge_range_id+1], sstruct_nodes)
                strands2 = domain_property.get_num_betasheet_strands( \
                        dom_ranges[merge_range_id+2], sstruct_nodes)

                if min(strands0, strands1, strands2) >= 4 and \
                   max(strands0, strands1, strands2) >= 5:
                    continue


                # create temporarily dom_ranges_test (from dom_ranges) by
                # merging 2 adjacent domain subsections
                dom_ranges_test = copy.deepcopy(dom_ranges)
                dom_ranges_test = data.merge_ranges(dom_ranges_test, \
                                    merge_range_id, merge_range_id+1)


                # create list merging_scores of merging candidates where at
                # least one of both merging candidates reaches the
                # minimum score
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                merging_scores, merging_score_relations, \
                    merging_score_bothabove, cache = \
                        merging_score.get_merging_scores( \
                            dom_ranges_test, helix_ranges, \
                            cache, params, sess)

                if merging_scores[merge_range_id] > -1:
                    len0 = data.range_len(dom_ranges_test[merge_range_id])
                    len1 = data.range_len(dom_ranges_test[merge_range_id+1])
                    len_prop = max(len0, len1)/min(len0, len1)

                    if len_prop <= 5:
                        merged_ranges = data.merge_ranges( \
                            dom_ranges_test, merge_range_id, merge_range_id+1)

                if len(merged_ranges) != dom_ranges_len:
                    break

            if len(merged_ranges) != dom_ranges_len:
                dom_ranges_len = len(merged_ranges)
                dom_ranges = merged_ranges
                remaining_iterations = dom_ranges_len
                break
            else:
                remaining_iterations -= 1


        remaining_iterations = 1000
        dom_ranges_len = len(dom_ranges)

        for iteration in range(0, len(dom_ranges)):
            ctl.p('iteration '+str(iteration))

            if remaining_iterations < 1:
                ctl.p('finished after '+str(iteration)+' iterations')
                break

            for merge_range_id in range(0, len(dom_ranges)-2):
                strands0 = domain_property.get_num_betasheet_strands( \
                        dom_ranges[merge_range_id], sstruct_nodes)
                strands1 = domain_property.get_num_betasheet_strands( \
                        dom_ranges[merge_range_id+1], sstruct_nodes)
                strands2 = domain_property.get_num_betasheet_strands( \
                        dom_ranges[merge_range_id+2], sstruct_nodes)

                if min(strands0, strands1, strands2) >= 4 and \
                   max(strands0, strands1, strands2) >= 5:
                    continue


                # create temporarily dom_ranges_test (from dom_ranges) by
                # merging 2 adjacent domain subsections
                dom_ranges_test = copy.deepcopy(dom_ranges)
                dom_ranges_test = data.merge_ranges(dom_ranges_test, \
                                    merge_range_id+1, merge_range_id+2)


                # create list merging_scores of merging candidates where at
                # least one of both merging candidates reaches the
                # minimum score
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                merging_scores, merging_score_relations, \
                    merging_score_bothabove, cache = \
                        merging_score.get_merging_scores( \
                            dom_ranges_test, helix_ranges, \
                            cache, params, sess)

                if merging_scores[merge_range_id] > -1:
                    len0 = data.range_len(dom_ranges_test[merge_range_id])
                    len1 = data.range_len(dom_ranges_test[merge_range_id+1])
                    len_prop = max(len0, len1)/min(len0, len1)

                    if len_prop <= 5:
                        merged_ranges = data.merge_ranges( \
                            dom_ranges_test, merge_range_id, merge_range_id+1)

                if len(merged_ranges) != dom_ranges_len:
                    break

            if len(merged_ranges) != len(dom_ranges):
                dom_ranges_len = len(merged_ranges)
                dom_ranges = merged_ranges
                remaining_iterations = dom_ranges_len
                break
            else:
                remaining_iterations -= 1


    # step 7:
    # check for incomplete domains and try to merge them
    # For each incomplete domain, perform a merge if at least one of both
    # merging candidates reaches a reduced the minimum score.
    # ===================================================================
    incomplete_dom_ranges, incomplete_dom_range_ids = \
            check_domain.check_domain_properties( \
                dom_ranges, sstruct_nodes, helix_ranges)

    merging_scores, merging_score_relations, merging_score_bothabove, \
        cache = \
            merging_score.get_merging_scores( \
                dom_ranges, helix_ranges, cache, params, sess, 0, 0.4)

    ranges_to_merge = []

    for incomplete_dom_range_id in incomplete_dom_range_ids:
        incomplete_dom_range = dom_ranges[incomplete_dom_range_id]

        max_score_merge_partner_score = -1
        max_score_merge_partner_range_id = -1

        if incomplete_dom_range_id > 0:
            if merging_scores[incomplete_dom_range_id-1] > \
               max_score_merge_partner_score:
                max_score_merge_partner_score = \
                        merging_scores[incomplete_dom_range_id-1]
                max_score_merge_partner_range_id = incomplete_dom_range_id-1

        if incomplete_dom_range_id < len(dom_ranges)-1:
            if merging_scores[incomplete_dom_range_id] > \
               max_score_merge_partner_score:
                max_score_merge_partner_score = \
                        merging_scores[incomplete_dom_range_id]
                max_score_merge_partner_range_id = incomplete_dom_range_id

        if max_score_merge_partner_score > -1:
            ranges_to_merge.append( (max_score_merge_partner_range_id, \
                                     max_score_merge_partner_range_id+1) )


    dom_ranges = data.merge_multiple_ranges(dom_ranges, ranges_to_merge)


    if options['verbous']:
        ctl.p('domain section ranges after merging:')
        ctl.p(dom_ranges)


    # transfer merging results to dssp_matr2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sepcand_sums = data.ranges_to_mask(dom_ranges, seq_len)
        # first res is 0

    for x in range(0, seq_len):
        for y in range(0, seq_len):
            if sepcand_sums[x] != 0 and sepcand_sums[y] != 0:
                dssp_matr2[x][y] = 1


    # plot matrices
    if options['show_plot']:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))

        ax[0][0].imshow(dssp_matr_betasheet_interstrand)
        ax[0][1].imshow(dssp_matr0)
        ax[1][0].imshow(dssp_matr_betasheet_interstrand2)
        ax[1][1].imshow(dssp_matr2) # with crosslinks

        plt.tight_layout()
        plt.show()


    # color dom_ranges and "linker ranges" between dom_ranges
    res_to_col = ''
    for i,s in enumerate(sepcand_sums):
        if s == 0:
            res_to_col += str(i+1)+','

    sess.run('color #1.1:'+res_to_col[:-1]+' red')
    visual.color_dom_sections(dom_ranges, 'gray', \
                              options['dom_coloring'], sess)
    if options['show_steps']:
        model_id += 1
        sess.run('combine #1 close false modelId #'+str(model_id))
        sess.run('rename #'+str(model_id)+' lowering')


    # step 8:
    # postprocessing / further versions
    # =================================

    # concatenate small domains with small alpha helices
    # --------------------------------------------------
    dom_ranges = specialdomains.merge_alphahelixdomains( \
                        dom_ranges, helix_ranges, sstruct_nodes, (20, 40))
    dom_ranges = specialdomains.merge_alphahelixdomains( \
                        dom_ranges, helix_ranges, sstruct_nodes, (20, 40))

    sess.run('color #1.1:'+res_to_col[:-1]+' red')
    visual.color_dom_sections(dom_ranges, 'gray', \
                              options['dom_coloring'], sess)
    if options['show_steps']:
        model_id += 1
        sess.run('combine #1 close false modelId #'+str(model_id))
        sess.run('rename #'+str(model_id)+' "small helices"')


    # version with domains containing part of linker
    # ----------------------------------------------

    # dom_boundaries: include half of linker
    dom_boundaries = [[] for i in range(0,len(dom_ranges))]
    dom_boundaries[0] = [1, seq_len]

    for i,d in enumerate(dom_ranges):

        if i < len(dom_ranges)-1:
            linker_half_pos, cache['linker_midpoints'] = \
                linker.get_linkermidpoint( \
                    [dom_ranges[i], dom_ranges[i+1]], \
                    cache['linker_midpoints'], sess)
            dom_boundaries[i+1] = [linker_half_pos+1, seq_len]
        else:
            linker_half_pos = seq_len

        dom_boundaries[i][1] = linker_half_pos


    # merging of domain sections at the N terminal end
    # ------------------------------------------------

    # merge potential signal sequence with adjacent range
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while specialdomains.has_separate_signal_sequence_dom(dom_ranges):
        dom_ranges = data.merge_ranges(dom_ranges, 0, 1)
        dom_boundaries = data.merge_ranges(dom_boundaries, 0, 1)


    # merge potential N terminal domain section containing a single betastrand
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if specialdomains.has_nterm_single_betastrand(dom_ranges, helix_ranges, \
                                                  sstruct_nodes):
        dom_ranges = data.merge_ranges(dom_ranges, 0, 1)
        dom_boundaries = data.merge_ranges(dom_boundaries, 0, 1)


    # mark first residue of each chain
    # --------------------------------
    for d in dom_boundaries:
        sess.run('color #1.1:'+str(d[0])+'-'+str(d[0])+' blue')

    if options['show_steps']:
        sess.run('tile')

    if options['export_cxs']:
        sess.run('save "'+f[:-4]+'.cxs'+'"')

    if options['verbous']:
        ctl.p(dom_ranges)
        ctl.p('domain boundaries:')
        ctl.p(dom_boundaries)


    # export fasta file with domain separations
    # -----------------------------------------
    f_out = open(f[:-4]+'_d.fa', 'w')
    fasta.fasta_separations(f_out, seq, dom_boundaries)

    return
