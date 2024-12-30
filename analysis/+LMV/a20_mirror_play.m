%% Movies showing clusters of neuronal responses tiling various task phases

movDir = LMV.Data.GetAnalysisDir('movies', 'mirror_play');

%% Specify example unit

expId = 2;

switch expId
    case 1
        expClusId = 520200153;
        srFile = "sr_NP52_B2_2023-03-31_16-51-11_final-v2.mat";
        datPath = "D:\sorting_workspace\NP52_B2_g0_imec0\temp_wh.dat";
        senInd = 5;
    case 2
        expClusId = 460100003;
        srFile = "sr_NP46_B1_2023-02-06_22-53-06_final-v2.mat";
        datPath = "D:\sorting_workspace\NP46_B1_g0_imec0\temp_wh.dat";
        senInd = 5;
end

%% Load data

% Load se
recId = NP.SE.GetID(expClusId);
srcTb = LMV.Data.FindSource(recId);
se = NP.SE.LoadSession(srcTb.path);

% Remove other units
isUnit = expClusId == NP.Unit.GetClusTb(se).clusId;
NP.Unit.RemoveUnits(se, ~isUnit);

% Preprocess se
se = LMV.SE.Transform(se, "enrich");

% Load sorting result object (for extracting spiking audio)
load(fullfile(LMV.Data.GetAnalysisDir, "misc", "data", srFile), "sr");
sr.mdat.Filename = datPath;

%% Sentence time alignment (M1)

% Find morph times to the best matched trials
NP.SE.SetMorphTimes(se, [], 'best');

% Morph se
seMorph = NP.SE.MorphSession(se);

%% Split stim and prod phases to separate se objects

% Reslice se at stim and prod cue onsets
[tt, tv] = seMorph.GetTable('taskTime', 'taskValue');
rt = seMorph.GetReferenceTime('taskTime');
tCue = sort([tt.cue1On+rt; tt.cue3On+rt]);
seRS = seMorph.Duplicate;
seRS.SliceSession(tCue, 'absolute');

% Unify event times
seRS.SetColumn('taskTime', 'trialOn', zeros(seRS.numEpochs, 1));
tt = seRS.GetTable("taskTime");
seRS.SetColumn('taskTime', 'speechMatchOn', max(tt.stimMatchOn, tt.prodMatchOn));
seRS.SetColumn('taskTime', 'speechMatchOff', max(tt.stimMatchOff, tt.prodMatchOff));

% Separate epochs
se2 = seRS.Split({1:2:seRS.numEpochs, 2:2:seRS.numEpochs});
phaseNames = ["stim", "prod"];
for i = 1 : numel(se2)
    tv.phase(:) = phaseNames(i);
    se2(i).SetTable("taskValue", tv, 'eventValues');
end

%% Compute PETHs and make sentence tables

[ce2, senTb2] = arrayfun(@(x) LMV.SE.ComputeSentencePETH(x, 'BoundaryEvents', {'trialOn', 'speechMatchOff'}), se2, 'Uni', false);

% % Find peak spike rate for each phase and each unit
% for i = 1 : numel(ce2)
%     rMax = max(arrayfun(@(x) x.clusTb.peakSpkRate, senTb2{i}.ce));
%     ce2{i}.clusTb.peakSpkRate = rMax;
%     for j = 1 : height(senTb2{i})
%         senTb2{i}.ce(j).clusTb.peakSpkRate = rMax;
%     end
% end

% Align to matched speech onsets
for i = 1 : numel(ce2)
    ce2{i}.AlignTime("speechMatchOn", "taskTime");
    arrayfun(@(x) x.AlignTime("speechMatchOn", "taskTime"), senTb2{i}.se);
end

% Put template trials to a new se (sentences should be in the same order as those in ce)
seTemp2 = cell(size(senTb2));
for j = 1 : numel(seTemp2)
    senTb = senTb2{j};
    clear seTemp
    for i = height(senTb) : -1 : 1
        m = senTb.trialNum{i} == senTb.tempTrialNum(i);
        seTemp(i) = senTb.se(i).Split({m});
    end
    seTemp = Merge(seTemp);
    seTemp.userData = seTemp.userData(1);
    seTemp2{j} = seTemp;
end

% Custom order
if expClusId == 520200153
    for j = 1 : numel(senTb2)
        senTb2{j}.se(senInd).SortEpochs([1:4 6 7 5]);
    end
end

%% 

for i = 1 : numel(ce2)
    pn = phaseNames(i);
    
    ce2{i}.userData.phaseName = pn;
    ce2{i}.userData.expClusTb = ce2{i}.clusTb;
    
    for j = 1 : height(senTb2{i})
        senTb2{i}.phase(j) = pn;
        senTb2{i}.se(j).userData.expClusTb = ce2{i}.clusTb;
    end
    
    seTemp2{i}.userData.phaseName = pn;
    seTemp2{i}.userData.expClusTb = ce2{i}.clusTb;
end

[ceStim, ceProd] = ce2{:};
[senTbStim, senTbProd] = senTb2{:};
[seTempStim, seTempProd] = seTemp2{:};
suffixes = ["Stim", "Prod"];

%% Make movies for stim

sf = "Stim";

if ~exist("mp"+sf, 'var')
    mpStim = MPlotter();
end
LMV.Linker.Mov.SetupPlotter(mpStim, sf);
            
mpStim.timeLimits = [-0.3 2];

%% 

LMV.Linker.Mov.MakeVideos(movDir, mpStim, seTempStim, senInd);

%% Make speech audio

LMV.Linker.Mov.MakeSpeechAudio(movDir, seTempStim, senInd);

%% Make spike audio

LMV.Linker.Mov.MakeSpikeAudio(movDir, seTempStim, senInd, se, sr);


%% Make movies for prod

sf = "Prod";

if ~exist("mp"+sf, 'var')
    mpProd = MPlotter();
end
LMV.Linker.Mov.SetupPlotter(mpProd, sf);
            
mpProd.timeLimits = [-0.3 3];

%% 

LMV.Linker.Mov.MakeVideos(movDir, mpProd, seTempProd, senInd);

%% Make speech audio

LMV.Linker.Mov.MakeSpeechAudio(movDir, seTempProd, senInd);

%% Make spike audio

LMV.Linker.Mov.MakeSpikeAudio(movDir, seTempProd, senInd, se, sr);



