# Activate the aligner environment
source /opt/anaconda3/bin/activate /userdata/qgreicius/conda_envs/mfa
# Set the lib path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/userdata/qgreicius/conda_envs/mfa/lib/

# RUN VALIDATION...deal with any errors if they arise before you move onto the next bit.
/userdata/qgreicius/conda_envs/mfa/bin/mfa validate /userdata/qgreicius/Neuropixels/preproc/NP59_B1/phones spanish_mfa spanish_mfa

mfa configure retry_beam 100
mfa configure --always_overwrite
# RUN ALIGNMENT
mfa align /userdata/qgreicius/Neuropixels/preproc/NP59_B1/phones spanish_mfa spanish_mfa /userdata/qgreicius/Neuropixels/preproc/NP59_B1/phones/results --clean