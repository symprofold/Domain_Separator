import ctl
import chimerax_api
import data
import ses


def get_contact_separation_residues(range_, ranges_potential_contact, sess):
    '''
    Calculate contacts for separation at all separation residues in the list
    range_ result.
    '''
    model_id = (1, 1)
    contact_resids = []
    contact_resids_quan = []

    resids = data.ranges_to_resid_list([range_])


    # performance optimization: calculate contacts for last separation point
    # in range_
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    resid = resids[-1]

    resids0_str = ses.ranges_to_str( \
                    [[ranges_potential_contact[0][0], resid-1]])
    resids1_str = ses.ranges_to_str( \
                    [[resid+1, ranges_potential_contact[1][1]]] )

    sess.run('close #'+str(model_id[0])+'.'+str(model_id[1])+'.1')

    contacts = sess.run('contact #'+ \
            str(model_id[0])+'.'+str(model_id[1])+':'+resids0_str+ \
            ' restrict #'+ \
            str(model_id[0])+'.'+str(model_id[1])+':'+resids1_str)

    sess.run('close #'+str(model_id[0])+'.'+str(model_id[1])+'.1')

    if len(contacts) == 0:
        zero_at_end = True
    else:
        zero_at_end = False


    for i,resid in enumerate(resids):

        # performance optimization: lowering of sampling rate
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if i > 5 and len(contact_resids) <= 1 and i%2 != 0 and \
            len(resids) > 10:
            continue
        elif i > 10 and len(contact_resids) <= 2 and i%4 != 0 and \
             len(resids) > 20:
            continue
        elif i > 40 and len(contact_resids) <= 4 and i%8 != 0:
            continue
        elif i > 80 and len(contact_resids) <= 4 and i%16 != 0:
            continue
        elif i > 160 and len(contact_resids) <= 4 and i%32 != 0:
            continue


        resids0_str = ses.ranges_to_str( \
                [[ranges_potential_contact[0][0], resid-1]])
        resids1_str = ses.ranges_to_str( \
                [[resid+1, ranges_potential_contact[1][1]]])

        sess.run('close #'+str(model_id[0])+'.'+str(model_id[1])+'.1')

        contacts = sess.run('contact #'+ \
                str(model_id[0])+'.'+str(model_id[1])+':'+resids0_str+ \
                ' restrict #'+ \
                str(model_id[0])+'.'+str(model_id[1])+':'+ resids1_str)

        sess.run('close #'+str(model_id[0])+'.'+str(model_id[1])+'.1')


        # performance optimization: return zero contact residues if first and
        # last last separation point in range_ has no contacts
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if i == 0 and len(contacts) == 0 and zero_at_end:
            return [], []

    return contact_resids, contact_resids_quan


def calc_linkermidpoint(range_, ranges_potential_contact, sess):
    '''
    Determine midpoint of free linker or free loop.
    '''
    contact_resids, contact_resids_quan = get_contact_separation_residues( \
            range_, ranges_potential_contact, sess)
    
    max_diff = 1
    max_diff_range = range_

    for i in range(1, len(contact_resids)):
        diff = contact_resids[i]-contact_resids[i-1]
        if diff > max_diff:
            max_diff = diff
            max_diff_range = [contact_resids[i-1], contact_resids[i]]

    midpoint = round(((max_diff_range[1]+max_diff_range[0])/2)-0.49)
        # -0.49: to ensure desired rounding behaviour

    if max_diff == 1 and len(contact_resids_quan) > 1:
        contact_resids_quan_sort = \
                sorted(contact_resids_quan, key=lambda x: x[1])
        midpoint = contact_resids_quan_sort[0][0]

    if midpoint not in data.ranges_to_resid_list([range_]):
        ctl.e(range_)
        ctl.e(midpoint)
        ctl.error('calc_linkermidpoint: midpoint not in range')

    range_resids = data.ranges_to_resid_list([range_])

    if len(range_resids) > 1 and midpoint == range_resids[-1]:
        midpoint = range_resids[-2]
    elif len(range_resids) == 1 and midpoint == range_resids[0]:
        midpoint = range_resids[0]-1

    return midpoint


def get_linkermidpoint(ranges, linkermidpoints, sess):
    '''
    Lookup if midpoint of free linker or free loop was calculated before and
    calculate if neccessary.
    '''
    model_id = (1, 1)
    res_ranges_key = []

    for res_range in ranges:
        res_ranges_key += res_range

    key = model_id+tuple(res_ranges_key)

    if key in linkermidpoints:

        return linkermidpoints[key], linkermidpoints

    else:
        linkermidpoint = calc_linkermidpoint( \
                [ranges[0][1]+1, ranges[1][0]-1], ranges, sess)
        linkermidpoints[key] = linkermidpoint

        return linkermidpoint, linkermidpoints
