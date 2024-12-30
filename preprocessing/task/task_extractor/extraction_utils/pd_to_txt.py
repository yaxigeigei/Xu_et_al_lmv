import pandas as pd
def df_to_txt(df, fname):
    """
    Take a pandas dataframe
    
    Send it to a text file that is readable by audacity.
    """
    if fname is None:
        print('no filename specified, writing temporary output :sob:')
        fname = './tmp_lab.txt'
    with open(fname, 'w') as f:
        for on, off, label in zip(df['Onset time'], df['Offset time'], df['Type']):
            s = '\t'.join(['%.6f' %on, '%.6f' %off, label]) + '\n'
            f.write(s)
            
            
        