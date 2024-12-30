%% Parameters

% Recording info
subjId = "NP112";
blockId = "B1";
recId = subjId + '_' + blockId;

% KS defaults
ksArgs = struct;
ksArgs.ConfigFunc = @MKilosort2.Config384;
ksArgs.ChannelMapFile = []; % extract from ap.meta file by default

% Special channel maps
switch char(recId)
    case 'NP54_B1', ksArgs.ChannelMapFile = char(recId + "_kilosortChanMap.mat"); % excluded a range of bad channels
    case 'NP_organoid_B2', ksArgs.ChannelMapFile = char("NP_organoid_kilosortChanMap.mat");
end

%% Copy and preprocess recording to userdata

catgt_copy;
ks_map; % generate full spike and lfp maps

%% Apply MTracer motion correction

mtracer_mc;
ks_map; % generate full spike and lfp maps

%% Spike sorting

switch char(recId)
    case 'NP04_B2', ksArgs.trange = [41, 931];
    case 'NP34_B2', ksArgs.trange = [180, 1223];
    case 'NP34_B3', ksArgs.trange = [100, Inf];
    case 'NP35_B1', ksArgs.trange = [0, 425];
    case 'NP36_B1', ksArgs.trange = [80, Inf];
    case 'NP37_B1', ksArgs.trange = [45, Inf];
    case 'NP37_B2', ksArgs.trange = [0, 405];
    case 'NP38_B5', ksArgs.trange = [14, 1055];
    case 'NP38_B6', ksArgs.trange = [0, 607];
    case 'NP40_B1', ksArgs.trange = [170, 774];
    case 'NP40_B2', ksArgs.trange = [130, 501];
    case 'NP41_B1', ksArgs.trange = [200, 730];
    case 'NP41_B2', ksArgs.trange = [50, 315];
    case 'NP42_B1', ksArgs.trange = [0, 750];
    case 'NP44_B2', ksArgs.trange = [80, 1077];
    case 'NP44_B3', ksArgs.trange = [82, 780];
    case 'NP45_B1', ksArgs.trange = [118, 685];
    case 'NP45_B2', ksArgs.trange = [96, 695];
    case 'NP46_B1', ksArgs.trange = [280, 762];
    case 'NP46_B2', ksArgs.trange = [170, 480]; %orig
    % case 'NP46_B2', ksArgs.trange = [215, 480]; %for shj
    case 'NP47_B1', ksArgs.trange = [175, 689];
    case 'NP47_B3', ksArgs.trange = [0, 650];
    case 'NP48_B1', ksArgs.trange = [165, 641];
    case 'NP49_B1', ksArgs.trange = [60, 1246];
    case 'NP50_B2', ksArgs.trange = [130, 360];
    case 'NP50_B3', ksArgs.trange = [80, 1174];
    case 'NP51_B1', ksArgs.trange = [90, 647];
    case 'NP52_B1', ksArgs.trange = [128,785];
    case 'NP52_B2', ksArgs.trange = [36, 810];
    case 'NP53_B1', ksArgs.trange = [102,1171];
    case 'NP54_B1', ksArgs.trange = [250, 850];
    case 'NP55_B1', ksArgs.trange = [130, 1142];
    case 'NP55_B2', ksArgs.trange = [100, 752];
    case 'NP55_B3', ksArgs.trange = [0, 85];
    case 'NP56_B1', ksArgs.trange = [160, 1576];
    case 'NP58_B2', ksArgs.trange = [345, 1180];
    case 'NP58_B4', ksArgs.trange = [59, 140];
    case 'NP59_B1', ksArgs.trange = [175, 1271];
    case 'NP59_B2', ksArgs.trange = [137, 673];
    case 'NP60_B1', ksArgs.trange = [175, 990];
    case 'NP61_B1', ksArgs.trange = [0, 790];
    case 'NP62_B1', ksArgs.trange = [317, 805];
    case 'NP62_B2', ksArgs.trange = [90, 817];
    case 'NP64_B1', ksArgs.trange = [113, 1249];
    case 'NP65_B2', ksArgs.trange = [0, 685]; %for task
    % case 'NP65_B2', ksArgs.trange = [100,400]; %for SHJ lab
    case 'NP66_B1', ksArgs.trange = [288, 1347];
    case 'NP66_B2', ksArgs.trange = [240, 1609];
    case 'NP67_B1', ksArgs.trange = [175, 2600];
    case 'NP68_B1', ksArgs.trange = [65, 415];
    case 'NP69_B1', ksArgs.trange = [173, 1277];
    case 'NP69_B2', ksArgs.trange = [87, 1029];
    case 'NP70_B2', ksArgs.trange = [0, 684];
    case 'NP70_B3', ksArgs.trange = [86, 969];
    case 'NP71_B2', ksArgs.trange = [0, 1431];
    case 'NP71_B3', ksArgs.trange = [77, 263];
    case 'NP72_B1', ksArgs.trange = [257, 1543];
    case 'NP72_B2', ksArgs.trange = [79, 1103];
    case 'NP73_B1', ksArgs.trange = [196, 846];
    case 'NP73_B2', ksArgs.trange = [143, 922];
    case 'NP74_B1', ksArgs.trange = [178, 905];
    % case 'NP74_B1', ksArgs.trange = [278, 906]; %for SHJ
    case 'NP75_B1', ksArgs.trange = [101, 875];
    case 'NP75_B2', ksArgs.trange = [169, 425];
    case 'NP76_B1', ksArgs.trange = [130, 1587];
%     case 'NP76_B1', ksArgs.trange = [530, 1587]; %for SHJ
    case 'NP77_B1', ksArgs.trange = [64, 234];
    case 'NP77_B2', ksArgs.trange = [0, 1062];
    case 'NP78_B1', ksArgs.trange = [523, 2156];
    case 'NP78_B2', ksArgs.trange = [65, 435];
    case 'NP79_B1', ksArgs.trange = [108, 630];
    case 'NP79_B2', ksArgs.trange = [250, 1988];
    case 'NP80_B1', ksArgs.trange = [148, 1577];
    case 'NP81_B1', ksArgs.trange = [120, 788];
    case 'NP81_B2', ksArgs.trange = [117, 1285];
    case 'NP82_B1', ksArgs.trange = [229, 639];
    case 'NP82_B3', ksArgs.trange = [303, 900];
    case 'NP82_B5', ksArgs.trange = [140, 712];
    case 'NP84_B1', ksArgs.trange = [113, 2500];
    case 'NP85_B1', ksArgs.trange = [120, 1000];
    case 'NP85_B2', ksArgs.trange = [10, 377];
    case 'NP85_B3', ksArgs.trange = [69, 1063];
    case 'NP86_B1', ksArgs.trange = [207, 900];
    case 'NP87_B1', ksArgs.trange = [113, 1730];
    case 'NP88_B1', ksArgs.trange = [0, 531];
    case 'NP88_B2', ksArgs.trange = [71, 1045];
    case 'NP89_B1', ksArgs.trange = [188, 1511];
    case 'NP89_B2', ksArgs.trange = [32, 332];
    case 'NP90_B1', ksArgs.trange = [256, 2025];
    % case 'NP90_B1', ksArgs.trange = [456, 2025]; %reduced range for SHJ
    case 'NP90_B2', ksArgs.trange = [205, 852];
    case 'NP90_B3', ksArgs.trange = [15, 213];
    case 'NP90_B4', ksArgs.trange = [0, 411];
    case 'NP91_B1', ksArgs.trange = [238, 1560];
    case 'NP92_B1', ksArgs.trange = [660, 1375];
    case 'NP92_B2', ksArgs.trange = [131, 1721];
    case 'NP93_B1', ksArgs.trange = [213, 1292];
    case 'NP94_B1', ksArgs.trange = [120, 2079];
    case 'NP94_B3', ksArgs.trange = [0, 271];
    case 'NP95_B1', ksArgs.trange = [179, 1172];
    case 'NP96_B1', ksArgs.trange = [164, 1047];
    case 'NP97_B1', ksArgs.trange = [169, 1887];
    case 'NP98_B1', ksArgs.trange = [202, 1766];
    case 'NP99_B1', ksArgs.trange = [140, 1317];
    case 'NP100_B1', ksArgs.trange = [0, 2147];
    case 'NP100_B2', ksArgs.trange = [80, 400];
    case 'NP101_B1', ksArgs.trange = [158, 936];
    case 'NP101_B2', ksArgs.trange = [186, 500];
    case 'NP101_B3', ksArgs.trange = [102, 584];
    case 'NP102_B2', ksArgs.trange = [101, 1340];
    case 'NP102_B3', ksArgs.trange = [115, 386];
    case 'NP104_B2', ksArgs.trange = [86, 1415];
    case 'NP105_B1', ksArgs.trange = [295, 1359];
    case 'NP105_B2', ksArgs.trange = [233, 931];
    case 'NP106_B1', ksArgs.trange = [170, 2210];
    case 'NP108_B1', ksArgs.trange = [127, 1353];
    case 'NP110_B1', ksArgs.trange = [249, 407];
    case 'NP110_B2', ksArgs.trange = [130, 881];
    case 'NP111_B1', ksArgs.trange = [364, 2297];
    case 'NP112_B1', ksArgs.trange = [0, 1473];
    case 'NP113_B1', ksArgs.trange = [118, 2074];
    case 'NP114_B1', ksArgs.trange = [225, 2603];
    case 'NP115_B2', ksArgs.trange = [45, 1039];
    case 'NP116_B1', ksArgs.trange = [71, 934]; % first half of the recording
    case 'NP116_B2', ksArgs.trange = [1269, 1870]; % second half of the recording
    case 'NP117_B1', ksArgs.trange = [186, 2166];
    case 'NP118_B2', ksArgs.trange = [89, 1347];
    case 'NP119_B1', ksArgs.trange = [51, 2568];
    case 'NP120_B1', ksArgs.trange = [406, 3671];
    case 'NP121_B1', ksArgs.trange = [271, 756];
    case 'NP121_B2', ksArgs.trange = [192, Inf];
    case 'NP121_B3', ksArgs.trange = [0, 382];
    case 'NP121_B4', ksArgs.trange = [77,715];
    case 'NP121_B5', ksArgs.trange = [0, 630];
    case 'NP122_B1', ksArgs.trange = [310, 2598];
    case 'NP_organoid_B2', ksArgs.trange = [153, 870];
end

switch char(recId)
    % Turn NP76_B1 automatic mc off for imec0 only. 
    case {'NP44_B3', 'NP49_B1', 'NP50_B3', 'NP66_B2', 'NP76_B1', 'NP85_B3','NP93_B1'}
        ksArgs.nblocks = 0;                             
end
disp(ksArgs);

ks_sort;

%% Plot spike maps

plot_maps;

figDir = fileparts(fileparts(rezTb.folder{1}));
timestamp = datestr(now, 'yyyymmdd_HHMMSS'); % add timestamp to the filename
saveas(f, fullfile(figDir, sprintf("%s_spike_maps_%s.png", recId, timestamp)));

