%% Paths

addpath("/userdata/dxu/project_np/code/np_preproc/task/intensity");
addpath("/userdata/dxu/project_np/code/babble");
addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab/npy-matlab"));
addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab/ManyFunMatlab"));

% Folder path of the chunked wav files
% e.g. /userdata/dxu/project_np/preproc/NP41_B1/audio_files/clips
wavDir = "{wav_dir}";
wavSearch = MBrowse.Dir2Table(fullfile(wavDir, '*.wav'));

% Output folder of extracted pitch data
% e.g. /userdata/dxu/project_np/preproc/NP41_B1/intensity
outputDir = "{output_dir}";

%% Compute intensity features

fprintf("\nExtracting intensity features using find_peakRate from files in:\n%s\n\n", wavDir);
fprintf("Output files will be saved at:\n%s\n\n", outputDir);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for i = 1 : height(wavSearch)
    % Read audio waveform
    wavFile = fullfile(wavSearch.folder{i}, wavSearch.name{i});
    [w, fs] = audioread(wavFile);
    
    % Compute intensity features
    [env, peakRate, peakEnv] = find_peakRate(w, fs);
    env = env(:);
    peakRate = peakRate(:);
    peakEnv = peakEnv(:);
    
    % Save as mat and npy files
    matFile = fullfile(outputDir, replace(wavSearch.name{i}, '.wav', '.mat'));
    save(matFile, 'env', 'peakRate', 'peakEnv');
    
    npyFile = fullfile(outputDir, replace(wavSearch.name{i}, '.wav', '.npy'));
    writeNPY([env, peakEnv, peakRate], npyFile);
    
    fprintf("Processed %i/%i files\n", i, height(wavSearch));
end

disp("Done!")

