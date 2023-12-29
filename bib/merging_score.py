import data
import merging_area
import ses


def get_merging_score(dom_ranges, range_ids0, range_ids1, helix_ranges, \
                      cache, params, sess, \
                      lowering_fact=1):
    '''
    Calculate merging score for merging two domain subsections.
    Each domain subsections is defined as a list of range ids.
    '''
    range0_resids = data.ranges_to_resid_list( \
                            [dom_ranges[i] for i in range_ids0])
    range1_resids = data.ranges_to_resid_list( \
                            [dom_ranges[i] for i in range_ids1])
    helix_resids = data.ranges_to_resid_list(helix_ranges)

    fract, area_diff0, area0, area1, \
            cache['surface_area_dict'], cache['linker_midpoints'] = \
            ses.merging_area_diff( \
                (1,1), dom_ranges, range_ids0, range_ids1, \
                cache['linker_midpoints'], cache['surface_area_dict'], sess)

    range0_len = len(range0_resids)
    range1_len = len(range1_resids)

    merging_area0 = merging_area.get_min_contact_area( \
                        area_diff0/fract, params, range0_len)
    merging_area1 = merging_area.get_min_contact_area( \
                        area_diff0/fract, params, range1_len)

    area_diff = area_diff0/2
    score0 = area_diff/range0_len
    score1 = area_diff/range1_len

    range0_helix_resids = set(range0_resids) & set(helix_resids)
    range1_helix_resids = set(range1_resids) & set(helix_resids)

    sstruct_fract, sstruct_area_diff, cache['surface_area_dict'] = \
            ses.merging_contact_diff( \
                sess, (1,1), \
                [list(range0_helix_resids), list(range1_helix_resids)], \
                cache['surface_area_dict'])

    overall_score = score0*score1
    overall_score_bothabove = score0*score1
    overall_score_relation = min(score0, score1)/max(score0, score1)
    
    if area_diff >= merging_area0*lowering_fact or \
       area_diff >= merging_area1*lowering_fact or \
       (sstruct_fract/2 >= params['mincontactarea_alphahelix_fract'] and \
       sstruct_area_diff/2 >= params['mincontactarea_alphahelix']):
        pass
    else:
        overall_score = -1

    if area_diff >= merging_area0*lowering_fact and \
       area_diff >= merging_area1*lowering_fact:
        pass
    else:
        overall_score_bothabove = -1

    return overall_score, overall_score_relation, overall_score_bothabove, \
           cache


def get_merging_scores(dom_ranges, helix_ranges, \
                       cache, params, sess, gap_ranges_num=0, \
                       lowering_fact=1):
    '''
    Get merging scores for all domain section ranges.
    '''
    merging_scores = [-1 for i in range(0, len(dom_ranges)-1)]
    overall_score_relations = [-1 for i in range(0, len(dom_ranges)-1)]
    merging_scores_bothabove = [-1 for i in range(0, len(dom_ranges)-1)]

    for i in range(0, len(dom_ranges)-1-gap_ranges_num):

        if gap_ranges_num > 0:
            len0 = dom_ranges[i][1]-dom_ranges[i][0]+1
            gap_len = dom_ranges[i+gap_ranges_num][1]- \
                      dom_ranges[i+gap_ranges_num][0]+1
            len1 = dom_ranges[i+gap_ranges_num+1][1]- \
                   dom_ranges[i+gap_ranges_num+1][0]+1

            if gap_len > min(len0, len1):
                continue

        merging_score, overall_score_relation, \
            merging_score_bothabove, cache = \
            get_merging_score(dom_ranges, [i], [i+gap_ranges_num+1], \
                              helix_ranges, cache, \
                              params, sess, lowering_fact)

        merging_scores[i] = merging_score
        overall_score_relations[i] = overall_score_relation
        merging_scores_bothabove[i] = merging_score_bothabove

    return merging_scores, overall_score_relations, \
           merging_scores_bothabove, cache


def get_highest_score(merging_scores):
    '''
    Get highest merging score.
    '''
    range_id_max = -1
    score_max = -1

    for range_id,score in enumerate(merging_scores):
        if score > score_max:
            score_max = score
            range_id_max = range_id

    return range_id_max, score_max
