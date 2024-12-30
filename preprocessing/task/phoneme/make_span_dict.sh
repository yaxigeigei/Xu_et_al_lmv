# Activate the aligner environment
source /opt/anaconda3/bin/activate /home/smetzger/.conda/envs/aligner
# Set the lib path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/smetzger/.conda/envs/aligner/lib/ 

# RUN VALIDATION...deal with any errors if they arise before you move onto the next bit.
/home/smetzger/.conda/envs/aligner/bin/mfa g2p spanish_g2p /userdata/smetzger/prod_feat_extract/phoneme_extract/vocab.txt /userdata/smetzger/prod_feat_extract/phoneme_extract/spanish_dict.txt