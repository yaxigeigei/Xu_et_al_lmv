# Activate the aligner environment
source /opt/anaconda3/bin/activate /home/smetzger/.conda/envs/aligner
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/smetzger/.conda/envs/aligner/lib/ 

# RUN VALIDATION...deal with any errors if they arise before you move onto the next bit.
/home/smetzger/.conda/envs/aligner/bin/mfa validate /userdata/smetzger/neuropixels/phones/NP11_B4/ /userdata/smetzger/prod_feat_extract/phoneme_extract/custom_lexicon.txt english

mfa configure retry_beam 100
mfa configure --always_overwrite
# RUN ALIGNMENT
mfa align /userdata/smetzger/neuropixels/phones/NP11_B4/ /userdata/smetzger/prod_feat_extract/phoneme_extract/custom_lexicon.txt english /userdata/smetzger/neuropixels/phones/NP11_B4/results/