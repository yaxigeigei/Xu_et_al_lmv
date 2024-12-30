%% Single-unit recording

% Unit quality control
LMV.a2_unit_qc;

% Physiological properties of units
LMV.a2_unit_wf;
LMV.a2_unit_isi;

%% LMV task

% Transform data
LMV.a3_morph_se;
LMV.a3_epoch_se;

% LMV phoneme statistics
LMV.a5_make_corpus_se;
LMV.a5_phone_stats;
LMV.a5_aai_stats;

% Task performance
LMV.a6_task_perf;

%% Task phase responsiveness

% Task phase responsiveness
LMV.a8_extract_resp;
LMV.a8_signrank;
LMV.a8_ttest;
LMV.a8_zeta;
% LMV.a8_benchmark_resp; % to be updated

% NMF clustering on phase responses
LMV.b7_nmfc;
LMV.b7_nmfc_fig;

% Task phase preferences
LMV.a8_venn;
% LMV.a8_phase_pref;
% LMV.b4_cla_phase; % population task phase decoding

%% Sentence selectivity

% Sentence selectivity tests
LMV.a9_sent_kw;
LMV.a9_kw_stats;
LMV.a9_kw_examples;

% Sentence decoding
LMV.a9_sent_cla;
LMV.a9_cla_perf; % LMV.a9_inspect;

% Sentence RSA
LMV.a9_sent_rsa;
LMV.a9_rsa_perf;

% % Sentence preference
% LMV.a14_sent_tuning;
% LMV.a14_sent_corr;

% % Population sentence decoding
% LMV.b4_make_ce;
% LMV.b4_cla_trial;
% LMV.b4_cla_perf;

%% Unit cache

LMV.a8_cache_unit; % make unit cache

%% Figure 1: Task phase tiling

% Recording sites on brain recon
LMV.a1_recon;

% PETH
LMV.a4_compute_peth;
LMV.a4_resp_heatmaps;
% LMV.a4_probe_peth;

% Example sentence rasters
LMV.a4_example_units;

% Functional embedding
LMV.a4_compute_umap;

% Functional-physiologial-anatomical distributions
LMV.a8_fpa;

% Movies
LMV.a7_m1_movie; % movie for M1 sentence responses
% LMV.a7_m2_movie; % movie for M2 session responses

%% Figures 2-3: Linker clusters

% Prepare input
LMV.a12_morph_se;       % morph prod to stim
LMV.a12_morph_qc;       % quality control of time morphing
LMV.a12_compute_peth;   % compute sentence PETHs
LMV.a12_cache_seq;      % cache peri-phone sequences by recordings (for plotting)

% Model linker transformation
LMV.a12_fit_lm;

% Process linking results
LMV.a12_hc;             % create clusTb with hierarchical clustering results
LMV.a12_m2_peth;        % add m2 PETHs to clusTb
LMV.a12_cache_combined; % cache combined results by recordings
LMV.a12_find_pos;       % add posTb to clusTb
LMV.a12_positions;      % extract linked positions

% Overall responses
LMV.a12_lib_profile;
% LMV.a12_lib_raster;
LMV.a12_lib_peth;

% Examples
LMV.a12_examples;

% Bridge
% LMV.a19_resp_heatmaps;
% LMV.a19_corr_delay;

% Movies
LMV.a20_mirror_play;
LMV.a20_bridge_play;
LMV.a20_mirror_tone;

%% Figure 4: Linker origins

LMV.a12_spatial;
LMV.a12_waveform;

%% Figure 5: Population activity

% Prepare response data
LMV.a16_compute_sen_resp;
LMV.a16_export_sen_resp;
LMV.a16_compute_trial_resp;

% Compute SCA and PCA
% ...\babble\lmv\sca\compute_sca.py

% Plot latent
LMV.a16_plot_lib;
LMV.a16_compute_metrics;
LMV.a16_latent;

% Loadings of mirror and bridge cells
LMV.a16_loading;

% Sentence decoding
LMV.a21_compute_trial_resp;
LMV.a21_resp_cla; % LMV.a21_resp_cla_jobs;
% LMV.a21_sca_cla;
LMV.a21_perf;

%% Figure 6: Speech encoding

% Prepare data
LMV.a10_make_ce; % feature and response timeseries
LMV.a13_cache_seq; % sequence aligned rasters/PETHs and features
LMV.a13_seq_lib;

% Features
LMV.a10_features;

% Fit sliding time RFs
LMV.a10_fit_st; % LMV.a10_fit_st_jobs;

% Examine time-locked feature encoding
LMV.a10_frac_seg;
LMV.a10_lib;
LMV.a10_top_rf;
LMV.a10_pop_rf;
LMV.a12_linker_rf;

% Sprectral example
recId = 'NP45_B2';
target = "stim";
clusId = 450200052;
LMV.a13_example_spectral;

% AKT example
recId = 'NP54_B1';
target = "prod";
clusId = 540100167;
LMV.a13_example_akt;

% Phone production example
recId = 'NP41_B1';
target = "prod";
clusId = 410100245;
LMV.a13_example_phone;

% Phone perception example
recId = 'NP44_B3';
target = "stim";
clusId = 440300585;
LMV.a13_example_phone;

% Sequence modulation
LMV.a13_seq_mod;

return
%% Contrastive learning

LMV.a18_cache_seq;
LMV.a18_extract_feat;

%% Sequence decoding

LMV.a17_make_ce;

%% Speech feature decoding

LMV.b4_cla_phase;
LMV.b4_cla_perf;

LMV.b4_reg_speech;
LMV.b4_reg_perf;

%% Sparsity

% Quantify phasic vs non-phasic responses
LMV.a11_phasic;
LMV.a11_phasic_dist;

% Quantify firing sparsity
LMV.a11_compute_peth;
LMV.a11_compute_sparsity;
LMV.a11_examples;
LMV.a11_depth;

%% TRF

% Fit TRFs
LMV.a10_fit_speech; % a10_fit_speech_jobs;
LMV.a10_trf_r;
LMV.a10_trf_lib;

% Examples
LMV.a10_examples;

% Cluster TRFs
LMV.a10_hc;
LMV.a10_hc_xall;

%% Spike-triggered feature analysis

LMV.a15_make_ce;
LMV.a15_extract_feat;
LMV.a15_feat_dim_re;

%% RSA

LMV.b8_make_ce;
LMV.b8_rsa_ctrl;
LMV.b8_rsa_evt;
LMV.b8_rsa_ts;
LMV.b8_rsa_select;

%% 

% 2024/04/25 Simon's Global keynote talk
LMV.c1_peth_heatmap;
LMV.a7_m2_movie;
LMV.c1_resp_rasters;
LMV.c1_mirror;

% 2024/06/26 AREADNE talk
LMV.c2_mel;


%% Utilities

cid = 4101e5+[249 405 138 463 245 224];

f = MPlot.Figure(111); clf
LMV.Overview.SessionFromCache(cid, 'DataSource', 'm2');
