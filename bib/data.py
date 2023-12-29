import ctl
import numpy as np


# data format conversions
# ~~~~~~~~~~~~~~~~~~~~~~~

def matr_to_ranges(matr):
    '''
    Conversion from matrix (first res is 0) to list of ranges
    (first res is 1).
    Each range is represented by a list with 2 entries:
    [first redidue, last residue].
    '''
    dom_ranges = []
    matr_sums = np.sum(matr, axis=1)

    in_range = False
    for i,r in enumerate(matr_sums):

        if r > 0 and not in_range:
            dom_ranges.append([i+1, -1])
            in_range = True
        elif r > 0 and in_range:
            continue
        elif r == 0 and in_range:
            dom_ranges[-1][1] = i
            in_range = False
        elif r == 0 and not in_range:
            continue

    if len(dom_ranges) > 0:
        if dom_ranges[-1][1] == -1:
            dom_ranges[-1][1] = len(matr_sums)

    return dom_ranges


def ranges_to_mask(ranges, seq_len):
    '''
    Conversion from ranges (start with 1) to mask matrix (start with 0).
    If a residue is part of a range, matrix element is 1, otherwise
    matrix element is 0.
    '''
    mask = np.zeros(seq_len)

    for range_ in ranges:
        for r in range(range_[0], range_[1]+1):
            mask[r-1] = 1

    return mask


def ranges_to_resid_list(ranges):
    '''
    Conversion from ranges (start with 1) to list or residue ids
    (start with 1).
    '''
    list_ = []

    for range_ in ranges:
        for r in range(range_[0], range_[1]+1):
            list_.append(r)

    return list_


def resids_to_ranges(resids):
    '''
    Conversion from resids (start with 1) to ranges (start with 1).
    '''
    mask = np.zeros((10000, 1))
        #start with 0

    for r in resids:
        mask[r-1][0] = 1

    ranges = matr_to_ranges(mask)

    return ranges


# range properties
# ~~~~~~~~~~~~~~~~
def range_len(range_):
    '''
    Get length of range.
    '''
    length = range_[1]-range_[0]+1

    return length


# filter ranges
# ~~~~~~~~~~~~~

def filter_range(range_, allowed_resids, max_ranges=1000):
    '''
    Filter a given range.
    '''
    residues = [r for r in range(range_[0], range_[1]+1)]
    residues_filtered = set(residues) & allowed_resids

    ranges_filtered = resids_to_ranges(residues_filtered)

    if len(ranges_filtered) == 0:

        return []

    if len(ranges_filtered) > max_ranges:
        ctl.e(range_)
        ctl.e(ranges_filtered)
        ctl.error('filter_range: too many separated ranges after filtering')

    return ranges_filtered


# merge ranges
# ~~~~~~~~~~~~

def merge_multiple_ranges(dom_ranges, ranges_to_merge):
    ''' Merge multiple adjacent ranges, input is list of range ids. '''

    merged_resids = set(ranges_to_resid_list(dom_ranges))

    for i,rangeid_pair in enumerate(ranges_to_merge):
        merged_resids |= set(ranges_to_resid_list( \
            merge_ranges(dom_ranges, rangeid_pair[0], rangeid_pair[1])))

    merged_ranges = resids_to_ranges(merged_resids)

    return merged_ranges


def merge_ranges(dom_ranges, range_id0_0, range_id1_0, gap=0):
    '''
    Merge ranges, input are range ids.
    Potential gaps between ranges are included.
    '''
    ranges_merged = []

    # checks
    if range_id0_0 == range_id1_0:
        ctl.error('data.merge_ranges: same range id')

    if abs(range_id1_0-range_id0_0) != 1+gap:
        ctl.error('data.merge_ranges: ranges not adjacent')

    range_id0 = min(range_id0_0, range_id1_0)
    range_id1 = max(range_id0_0, range_id1_0)

    for i,d in enumerate(dom_ranges):
        if i == range_id0:
            ranges_merged.append([dom_ranges[i][0], dom_ranges[i+1+gap][1]])
        elif i in range(range_id0, range_id1+1):
            continue
        else:
            ranges_merged.append(d)
    
    return ranges_merged
