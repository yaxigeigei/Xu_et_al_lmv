%% Generate a library of sequence plots

inputDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_seq");
outputDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_feat_cl");
sdSearch = MBrowse.Dir2Table(fullfile(inputDir, "*_seqData.mat"));

%% 

for recIdx = 2 : height(sdSearch)
    % Load cached sequence data
    recId = NP.SE.GetID(sdSearch.name{recIdx});
    disp(recId);
    s = load(fullfile(sdSearch.folder{recIdx}, sdSearch.name{recIdx}), "se", "feats", "seqData");
    seeds = s.feats;
    se = s.se;
    clusTb = NP.Unit.GetClusTb(se);
    
    % Load phoneme stRFs
    targets = LMV.TRF.targets;
    trfTb = LMV.TRF.LoadModels("phone_"+targets, recId);
    [~, trfInd] = MMath.SortLike(trfTb.clusId, clusTb.clusId);
    
    
    % Loop through targets
    for tIdx = 1 : numel(targets)
        % Check output
        target = targets(tIdx);
        wfPath = fullfile(outputDir, sprintf("%s_%s_waveforms.mat", target, recId));
        if exist(wfPath, "file")
            fprintf("Data for %s %s have already been exported.\n", recId, target);
            continue
        end
        
        % Extract responses and features from se
        C = s.seqData.(target){"phone"}{:,3};
        stRFs = trfTb.("phone_"+target)(trfInd);
        
        tMax = 1;
        switch target
            case 'stim'
                tWin = [-0.1 tMax];
                phaseName = 'stim';
            case 'prod'
                tWin = [-tMax 0.1];
                phaseName = 'prod';
            case 'feedback'
                tWin = [-0.1 tMax];
                phaseName = 'prod';
            otherwise
                error("'%s' is not a supported target name.", target);
        end
        
        for k = 1 : numel(C)
            seqTb = C{k}.seqTb;
            if isempty(seqTb)
                continue
            end
            seqTb(seqTb.nSample<3,:) = [];
            seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin);
            % seqTb = LMV.PE.AddFeat2SeqTb(seqTb, se, flip(-tWin), "mel", "mic");
            seqTb = LMV.PE.AddWaveform2SeqTb(seqTb, se, flip(-tWin));
            C{k}.seqTb = seqTb;
        end
        
        % 
        [featTb, respTb] = LMV.PE.MakeFlattenedTables(C, stRFs);
        time = featTb.tWaveform{1};
        speaker1 = featTb.speaker1;
        mic = featTb.mic;
        featTb(:,["tWaveform", "speaker1", "mic"]) = [];
        
        % Save mat files
        % save(fullfile(outputDir, sprintf("%s_%s.mat", target, recId)), 'featTb', 'respTb');
        save(wfPath, 'time', 'speaker1', 'mic');
        
        fDf = py.pandas.DataFrame(featTb(:,1:6));
        rDf = py.pandas.DataFrame(respTb);
        fDf.to_pickle(fullfile(outputDir, sprintf("%s_%s_seq.pkl", target, recId)));
        rDf.to_pickle(fullfile(outputDir, sprintf("%s_%s_resp.pkl", target, recId)));
        % py.numpy.save(fullfile(outputDir, sprintf("%s_%s_mel.npy", target, recId)), cat(3, featTb.mel{:}));
        
        % return
    end
end
