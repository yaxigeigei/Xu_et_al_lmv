import os
import copy
import pandas as pd
import numpy as np

def totsv(csvfp):
    """
    This is a utility to convert csvs from previous recordings
    into text files that can be combined in this code base.
    
    input: csv file
    output: none, but saves the csv as a txt file 
    """
    df = pd.read_csv(csvfp)    
    nudf = df[['Onset time', 'Offset time', 'Type']]
    new = nudf.to_csv(csvfp.replace('.csv', '.txt'), sep='\t', header=False, index=False)

def todf(txtlist):
    """
    Utility to transform a text file from AUDACITY 
    Into a dataframe with Onset time, Offset time, and Type [type = stim timing]
    """
    with open(txtlist) as f: 
        lines = f.readlines()
    lines = [l for l in lines if not l.startswith('\\')]
    txtlist= lines
    on, off, types = [], [],[]
    for t in txtlist: 
        
        tl = t.split('\t')
        try:
            on.append(float(tl[0]))
            off.append(float(tl[1]))
            types.append(str(tl[-1].replace('\n', '')))
        except Exception as e:
            print(e)
            print('couldnt deal with this line:', tl)
            print('it was not added to the combined dataframe and hence wont be processed')
        
    import pandas as pd
    
    return pd.DataFrame({
        'Onset time':on, 
        'Offset time':off, 
        'Type':types
    })

def label_combine(directory, task_list):
    
    """
    inputs: directory - the labels directory
    task_list: each task should have a corresponding directory in the labels directory
    
    The structure of the directory should be: 
    task
        - speech_x_timing.txt where x is the source: prod (pt), listen/stim (e.g. the task), experimenter, whatever. 
            These files will undergo feature extraction. 
        - cue_timing or photodiode_timing, or reading_timing etc. The first word will be taken as the source. 
            - if speech is not in the label name, then it will NOT get features extracted. 
    """
    combined_speech_labels = []
    combined_reading_labels = []
    combined_task_labels = []
    
    
    for task in task_list:
        task_dir = os.path.join(directory, task)
        if not os.path.exists(task_dir):
            continue
        
        task_all_files = os.listdir(task_dir)
        task_txt_files = [t for t in task_all_files if t.endswith('txt')]
        
        for file in task_txt_files:
            if 'speech' in file:
                if not file.endswith('timing.txt'):
                    print('ignoring file:', file)
                    continue
                else: 
                    source = file.split('_')[1]
                    
                timing_df = todf(os.path.join(directory, task, file))
                timing_df['Source'] = [source]*len(timing_df)
                timing_df['Task'] = [task]*len(timing_df)
                combined_speech_labels.append(timing_df)
            elif 'reading' in file:
                if not file.endswith('timing.txt'):
                    print('ignoring file:', file)
                    continue
                else: 
                    source = file.split('_')[1]
                    
                timing_df = todf(os.path.join(directory, task, file))
                timing_df['Source'] = [source]*len(timing_df)
                timing_df['Task'] = [task]*len(timing_df)
                combined_reading_labels.append(timing_df)  
            else: 
                source = file.split('_')[0]
                timing_df= todf(os.path.join(directory, task, file))
                timing_df['Source'] = [source]*len(timing_df)
                timing_df['Task'] = [task]*len(timing_df)
                combined_task_labels.append(timing_df)
                
    if len(combined_speech_labels) > 0:
        combined_speech_labels = pd.concat(combined_speech_labels)
        combined_speech_labels = combined_speech_labels.sort_values('Onset time')
        combined_speech_labels.to_csv(os.path.join(directory, 'combined_speech_labels.csv'))

    if len(combined_reading_labels) > 0:
        combined_reading_labels = pd.concat(combined_reading_labels)
        combined_reading_labels = combined_reading_labels.sort_values('Onset time')
        combined_reading_labels.to_csv(os.path.join(directory, 'combined_reading_labels.csv'))
        
    if len(combined_task_labels) > 0:
        combined_task_labels = pd.concat(combined_task_labels)
        combined_task_labels = combined_task_labels.sort_values('Onset time')
        combined_task_labels.to_csv(os.path.join(directory, 'combined_task_labels.csv'))
        
                    
    print('labels combined and saved...enjoy :D')
            
            
            
        
        
    
def rename_files(directory, ind_to_time):
    """
    Renames files so that they have the start and end time.
    
    input: directory: a directory of files that has text grids or acoustic feature arrays. 
    ind_to_time, a dictionary of what timepoints an index corresponds to. This will be set by the prep_alignment file
    during phoneme processing. It should go up in order with time. 
    
    output: none, but will rename the files in the directory to match 
    """
    
    warnflag = False # warn people the path didnt get properly remade. 
    for file in os.listdir(directory): 
        
        if file.endswith('TextGrid'):
            og_file = copy.deepcopy(file)
            file = file.replace('stim_', '')
            file = file.replace('.TextGrid', '')
            end = '.TextGrid'
            
        elif file.endswith('.mat'):
            og_file = copy.deepcopy(file)
            
            file = file.replace('.mat', '')
            end = '.mat'
            
        elif 'formants' in file:
            og_file = copy.deepcopy(file)
            file = file.split('.npy')[0].replace('formants_', '')
            end = '.npy'
            
        elif file.endswith('.npy'): 
            og_file = copy.deepcopy(file)
            file = file.replace('.npy', '')
            end = '.npy'
        else: 
            print('no way to rename file currently')
            continue
        
        try: 
            os.rename(os.path.join(directory, og_file), 
                      os.path.join(directory, '%.3f_%.3f' %(ind_to_time[int(file)][0], ind_to_time[int(file)][1]) + end))
        except Exception as e:
            if not warnflag:
                print(str(e))
                warnflag =True
                print('you may have already changed the files. if so ignore the error!')