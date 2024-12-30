# Activate the aligner environment
source /opt/anaconda3/bin/activate /home/smetzger/.conda/envs/aligner
# Set the lib path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/smetzger/.conda/envs/aligner/lib/ 

# RUN VALIDATION...deal with any errors if they arise before you move onto the next bit.
/home/smetzger/.conda/envs/aligner/bin/mfa validate {corpus_dir} {lexicon} {language}

mfa configure retry_beam 100
mfa configure --always_overwrite
# RUN ALIGNMENT
mfa train {corpus_dir} {lexicon} {output_dir}