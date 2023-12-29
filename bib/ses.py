import data
import linker


def ranges_to_str(ranges):
    '''
    Conversion from ranges (start with 1) to string.
    '''
    ranges_txt = []

    for range_ in ranges:
        ranges_txt.append(str(range_[0])+'-'+str(range_[1]))

    txt = ','.join(ranges_txt)

    return txt


def calc_surface_area(sess, model_id, res_ranges):
    ''' Calculate surface area caused only by given residue range. '''
    
    area = -1

    if len(res_ranges) == 0:

        return 0

    # remove possibly preexisting object
    sess.run('close #'+str(model_id[0])+'.'+str(model_id[1])+'.1')

    ranges_str = ranges_to_str(res_ranges)

    sess.run('surface #'+ \
            str(model_id[0])+'.'+str(model_id[1])+':'+ranges_str+ \
            ' enclose #'+ \
            str(model_id[0])+'.'+str(model_id[1])+':'+ranges_str)

    area = sess.run('measure area #'+ \
            str(model_id[0])+'.'+str(model_id[1])+'.1 includeMasked false')
    area = float(area)
    sess.run('close #'+str(model_id[0])+'.'+str(model_id[1])+'.1')

    return area


def surface_area(sess, model_id, res_ranges, surface_area_dict):
    '''
    Lookup if surface area was calculated before and calculate if neccessary.

    Args:
        res_ranges: list of residue ranges
            multiple ranges possible to get total area of disconnected
            structure parts.
    '''
    res_ranges_key = []

    for res_range in res_ranges:
        res_ranges_key += res_range

    key = model_id+tuple(res_ranges_key)

    if key in surface_area_dict:

        return surface_area_dict[key], surface_area_dict

    else:
        area = calc_surface_area(sess, model_id, res_ranges)
        surface_area_dict[key] = area

        return area, surface_area_dict


def get_range_plus_half_linker(dom_ranges, range_id, linkermidpoints, sess):
    '''
    Determine the residue range including half of the linker on n-terminal and
    c-terminal side.
    '''
    if range_id > 0:
        midpoint0, linkermidpoints = linker.get_linkermidpoint( \
            [dom_ranges[range_id-1], dom_ranges[range_id]], \
            linkermidpoints, sess)
        midpoint0 += 1
    else:
        midpoint0 = dom_ranges[range_id][0]

    if range_id < len(dom_ranges)-1:
        midpoint1, linkermidpoints = linker.get_linkermidpoint( \
            [dom_ranges[range_id], dom_ranges[range_id+1]], \
            linkermidpoints, sess)
    else:
        midpoint1 = dom_ranges[range_id][1]

    range_with_linker = [midpoint0, midpoint1]

    return range_with_linker, linkermidpoints


def merging_area_diff(model_id, dom_ranges, range_ids0, range_ids1, \
                      linkermidpoints, surface_area_dict, sess):
    '''
    Calculate area difference between separated domain sections and
    merged domain sections.
    '''
    area_diff_fraction = -1

    if range_ids0[0] > range_ids1[0]:
        ctl.error('ERROR: merging_area_diff: error in range order')

    # check if only 1 range per domain section
    if len(range_ids0) != 1 or len(range_ids1) != 1:
        ctl.error('merging_area_diff: more than 1 range per domain section')

    range_plus_half_linker0, linkermidpoints = \
        get_range_plus_half_linker(dom_ranges, range_ids0[0], linkermidpoints, sess)
    range_plus_half_linker1, linkermidpoints = \
        get_range_plus_half_linker(dom_ranges, range_ids1[0], linkermidpoints, sess)

    area0, surface_area_dict = surface_area( \
            sess, (1,1), [range_plus_half_linker0], surface_area_dict)
    area1, surface_area_dict = surface_area( \
            sess, (1,1), [range_plus_half_linker1], surface_area_dict)
    area_merged, surface_area_dict = surface_area( \
            sess, (1,1), [range_plus_half_linker0, range_plus_half_linker1], \
            surface_area_dict)

    area_diff = area0+area1-area_merged
    area_diff_fraction = area_diff/min(area0, area1)

    return area_diff_fraction, area_diff, area0, area1, \
           surface_area_dict, linkermidpoints


def merging_contact_diff(sess, model_id, ranges, surface_area_dict):
    '''
    Calculate area difference between two lists of ranges separated and
    both lists of ranges merged.
    '''
    area_diff_fraction = -1

    ranges0 = data.resids_to_ranges(ranges[0])
    ranges1 = data.resids_to_ranges(ranges[1])

    area0, surface_area_dict = surface_area( \
            sess, (1,1), ranges0, surface_area_dict)
    area1, surface_area_dict = surface_area( \
            sess, (1,1), ranges1, surface_area_dict)
    area_merged, surface_area_dict = surface_area( \
            sess, (1,1), ranges0+ranges1, surface_area_dict)

    area_diff = area0+area1-area_merged

    if area_diff > 0: # avoid division by 0
        area_diff_fraction = area_diff/min(area0, area1)
    else:
        area_diff_fraction = 0

    return area_diff_fraction, area_diff, surface_area_dict
