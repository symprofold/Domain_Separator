import ctl


def clean_path(path):
    """ Clean path string. """

    path_cleaned = []    
    path_arr = path.split('/')
    for p in path_arr:
        if p == '..':
            path_cleaned.pop()
        else:
            path_cleaned.append(p)

    path_cl = '/'.join(path_cleaned)

    return path_cl
