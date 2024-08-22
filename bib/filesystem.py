import ctl


def clean_path(path):
    ''' Clean path string. '''

    path_cleaned = []    
    path_arr = path.split('/')
    for p in path_arr:
        if p == '..':
            path_cleaned.pop()
        else:
            path_cleaned.append(p)

    path_cl = '/'.join(path_cleaned)

    return path_cl


def get_file(file, attempt=0):
    ''' Open file and create list with lines. '''

    attempt += 1
    ra = 5
 
    try:
        f = open(file, 'r')
        ra = f.readable()
        if ra == True:
            lin = f.readlines()
        else:
            ra = 5
        f.close()

    except IOError:
        ctl.d("ERROR: file:"+file+", attempt:"+str(attempt))
        ctl.d("retry")

        time.sleep(0.1*attempt)
        return get_file(file)

    if ra != True:
        ctl.d(ra) 

    return lin
