

def get_min_contact_area(area, params, dom_sect_len):
    ''' Get minimal contact area for domain section merging. '''

    a = params['min_contact_area_a']
    b = params['min_contact_area_b']

    min_area = area * ( a + b*dom_sect_len )

    return min_area
