import ctl
import domain_property


def check_domain_property(range_, sstruct_nodes, helices):
    '''
    Check reasons against domain property of a residue range.
    '''
    domprop = True

    betasheet_max4strands_num = domain_property.get_num_betasheet( \
                                    range_, sstruct_nodes, 4)
    betasheet_num = domain_property.get_num_betasheet(range_, sstruct_nodes)
    helix_num = domain_property.get_num_helices(range_, helices)

    if betasheet_max4strands_num == 1 and betasheet_num == 1 and \
       helix_num == 0:
        domprop = False

    if betasheet_max4strands_num == 0 and betasheet_num == 0 and \
       helix_num == 1 and range_[1]-range_[0]+1 <= 10:
        domprop = False

    return domprop


def check_domain_properties(ranges, sstruct_nodes, helices):
    '''
    Check reasons against domain property of each residue range in
    list 'ranges'.
    '''
    incomplete_res_ranges = []
    incomplete_range_ids = []

    for range_id,range_ in enumerate(ranges):

        if not check_domain_property(range_, sstruct_nodes, helices):
                incomplete_res_ranges.append(range_)
                incomplete_range_ids.append(range_id)

    return incomplete_res_ranges, incomplete_range_ids
