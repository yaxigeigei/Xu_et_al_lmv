# Activate the aligner environment
source /opt/anaconda3/bin/activate {aligner_env}
# Set the lib path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{aligner_env}/lib/

# RUN VALIDATION...deal with any errors if they arise before you move onto the next bit.
{aligner_env}/bin/mfa validate {corpus_dir} spanish_mfa spanish_mfa

mfa configure retry_beam 100
mfa configure --always_overwrite
# RUN ALIGNMENT
mfa align {corpus_dir} spanish_mfa spanish_mfa {output_dir} --clean