import os
this_dir = os.path.dirname(__file__)


def make_alignment_files(corpus_dir, lexicon_file, output_dir, language, artic_env):
    
    corpus_dir = str(corpus_dir)
    output_dir = str(output_dir)
    lexicon_file = str(lexicon_file)
    
    newlines = []
    with open(os.path.join(this_dir, 'run_aligner_template.sh')) as f: 
        lines = f.readlines()
        
    for l in lines: 
        newlines.append(process_line(l, corpus_dir, lexicon_file, output_dir, language, artic_env))
        
    with open(os.path.join(this_dir, 'run_aligner.sh'), 'w') as f: 
        f.writelines(newlines)


def make_spanish_alignment_files(corpus_dir, lexicon_file, output_dir, language, artic_env):
    
    corpus_dir = str(corpus_dir)
    output_dir = str(output_dir)
    lexicon_file = str(lexicon_file)
    
    newlines = []
    # with open(os.path.join(this_dir, 'span_aligner_template.sh')) as f: 
    with open(os.path.join(this_dir, 'span_aligner_template_qg.sh')) as f: 

        lines = f.readlines()
        
    for l in lines: 
        newlines.append(process_line(l, corpus_dir, lexicon_file, output_dir, language, artic_env))
        
    with open(os.path.join(this_dir, 'run_aligner_spanish.sh'), 'w') as f: 
        f.writelines(newlines)

        
def process_line(l, corpus_dir, lexicon, output_dir, language='english', artic_env=None): 
    l = l.replace('{corpus_dir}', corpus_dir)
    l = l.replace('{lexicon}', lexicon)
    l = l.replace('{output_dir}', output_dir)
    l = l.replace('{language}', language)
    l = l.replace('{aligner_env}', artic_env)
    return l
