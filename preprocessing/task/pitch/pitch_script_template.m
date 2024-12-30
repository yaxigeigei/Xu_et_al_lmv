%% Paths

addpath("/userdata/dxu/project_np/code/np_preproc");
addpath("/userdata/dxu/project_np/code/np_preproc/task/pitch");
addpath("/userdata/dxu/project_np/code/third_party_matlab/VoiceSauce");
addpath("/userdata/dxu/project_np/code/third_party_matlab/mPraat");
addpath("/userdata/dxu/project_np/code/third_party_matlab/npy-matlab");
addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab/ManyFunMatlab"));

% Folder path of the chunked wav files
% e.g. /userdata/dxu/project_np/preproc/NP41_B1/audio_files/clips
wavDir = "{wav_dir}";
wavSearch = MBrowse.Dir2Table(fullfile(wavDir, '*.wav'));

% Folder path of MFA TextGrid output
% e.g. /userdata/dxu/project_np/preproc/NP41_B1/phones/results/all
tgDir = "{phone_dir}";

% Output folder of extracted pitch data
% e.g. /userdata/dxu/project_np/preproc/NP41_B1/pitch/straight_output
outputDir = "{output_dir}";


%% Pitch estimation

fprintf("\nExtracting pitch using STRAIGHT from files in:\n%s\n\n", wavDir);
fprintf("Output files will be saved at:\n%s\n\n", outputDir);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Set up STRAIGHT parameters
ops.maxstrdur = 100;    % batch size in sec
ops.minstrF0 = 75;      % lower F0 bound in Hz
ops.maxstrF0 = 300;     % higher F0 bound in Hz

for i = 1 : height(wavSearch)
    % Read audio file
    wavFile = fullfile(wavSearch.folder{i}, wavSearch.name{i});
    [w, fs] = audioread(wavFile);
    
    % Run STRAIGHT
    [f0, v] = func_StraightPitch(double(w), fs, ops);
    t = (1 : numel(f0))' * 1e-3; % f0 is by default estimated at 1kHz
    
    % Refine pitch
    tgName = replace(wavSearch.name{i}, {'wav_', '.wav'}, {'stim_', '.TextGrid'});
    tgFile = fullfile(tgDir, tgName);
    if exist(tgFile, 'file')
        tg = NP.TGEvent(tgRead(tgFile));
        f0c = NP.Pitch.RefinePitch(f0, t, tg);
    else
        fprintf("Cannot find the TextGrid file to refine the pitch timeseries.\n");
        f0c = f0;
    end
    
    % Save original and refined pitch contour to an npy file
    npyName = strrep(wavSearch.name{i}, '.wav', '.npy');
    npyFile = fullfile(outputDir, npyName);
    writeNPY([f0 f0c], npyFile);
    
    % Downsample refined contour and save to a PitchTier file
    tq = t(1) : 5e-3 : t(end); % 5ms per sample
    f0cq = interp1(t, f0c, tq);
    isVal = ~isnan(f0cq); % PitchTier does not accept NaN
    
    pt = struct;
    pt.t = tq(isVal);
    pt.f = f0cq(isVal);
    
    ptName = strrep(wavSearch.name{i}, '.wav', '.PitchTier');
    ptFile = fullfile(outputDir, ptName);
    
    ptWrite(pt, ptFile);
    
    fprintf("Processed %i/%i files\n", i, height(wavSearch));
end

disp("Done!");

