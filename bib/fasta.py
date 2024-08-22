

def fasta_separations(export_file, seq, dom_ranges):
    '''
    Export to fasta file with domain separations.
    '''
    export_file.write('>sequence'+"\r\n")

    for d in dom_ranges:
        export_file.write(seq[d[0]-1:d[1]]+"\r\n\r\n")

    export_file.close()
    
    return
