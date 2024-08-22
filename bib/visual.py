import chimerax_api
import matplotlib.colors as mcolors


def color_dom_sections(dom_sections, color0, dom_coloring, sess): 
    ''' Color domain sections in visualization. '''

    colordict = {}

    for s in dom_sections:

        if dom_coloring:
            colordict, color_name = assign_color( \
                    colordict, str(s[0])+'-'+str(s[1]))
        else:
            color_name = color0
        
        sess.run('color #1.1:'+str(s[0])+'-'+str(s[1])+' '+str(color_name))

    return


def assign_color(colordict, sect_name):
    ''' Assign unique color to section. '''

    if sect_name in colordict:
        color_name = colordict[sect_name]
    else:
        color_names = dict(mcolors.TABLEAU_COLORS, **mcolors.CSS4_COLORS)
        colors_used = []

        for c in colordict:
            colors_used.append(colordict[c])

        for cn in color_names:

            if cn not in colors_used:
                color_name = cn
                colordict[sect_name] = cn
                break

    if color_name in color_names:
        color_name_return = color_names[color_name]
    else:
        color_name_return = 'black'

    return colordict, color_name_return
