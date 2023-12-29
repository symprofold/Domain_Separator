import data
import domain_property


def has_separate_signal_sequence_dom(dom_ranges):
    '''
    Check if list of domain ranges contains a separate N terminal signal
    sequence domain.
    '''
    if dom_ranges[0][1] <= 30:
        return True

    else:
        return False


def has_nterm_single_betastrand(dom_ranges, helices, sstruct_nodes):
    '''
    Check if list of domain ranges contains a N terminal single betastrand
    and maximal 1 alpha helices.
    '''
    strands_num = domain_property.get_num_betasheet_strands(dom_ranges[0], \
                                                            sstruct_nodes)
    helix_num = domain_property.get_num_helices(dom_ranges[0], helices)

    if strands_num == 1 and helix_num <= 1:

        return True

    else:

        return False


def check_alphahelixdomain_pair(range0, range1, helices, sstruct_nodes, \
                                max_range_lengths):
    '''
    Check if a pair of adjacent domain section ranges can be merged.
    The criterion is the presence of only alpha helices as
    secondary structure.

    Args:
        max_range_lengths: tuple with both maximum sequence lengths to merge
    '''
    helix0_num = domain_property.get_num_helices(range0, helices)
    betasheet0_num = domain_property.get_num_betasheet(range0, \
                                                      sstruct_nodes)

    helix1_num = domain_property.get_num_helices(range1, helices)
    betasheet1_num = domain_property.get_num_betasheet(range1, \
                                                      sstruct_nodes)

    if betasheet0_num == 0 and betasheet1_num == 0 and \
        helix0_num >= 1 and helix1_num >= 1 and \
        min(range0[1]-range0[0]+1, range1[1]-range1[0]+1) <= \
        min(max_range_lengths) and \
        max(range0[1]-range0[0]+1, range1[1]-range1[0]+1) <= \
        max(max_range_lengths):

        return True

    return False


def merge_alphahelixdomains(dom_ranges, helices, sstruct_nodes, \
                            max_range_lengths):
    '''
    Merge all possible pairs of adjacent domain section ranges that have only
    alpha helices as secondary structure.

    Args:
        max_range_lengths: tuple with both maximum sequence lengths to merge
    '''
    ranges_to_merge = []

    for range_id0 in range(0, len(dom_ranges)-1):

        if check_alphahelixdomain_pair( \
            dom_ranges[range_id0], dom_ranges[range_id0+1], \
            helices, sstruct_nodes, max_range_lengths):

            ranges_to_merge.append((range_id0, range_id0+1))

    merged_ranges = data.merge_multiple_ranges(dom_ranges, ranges_to_merge)

    return merged_ranges
