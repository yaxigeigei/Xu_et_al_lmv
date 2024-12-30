%% Initialize fearture table

tbPath = fullfile(LMV.Data.GetAnalysisDir, "encode", "feature_table_raw.xlsx");

artic = LMV.TRF.GetFeatureSet("artic3");
[articTb, articDict] = NP.Artic.GetFeatureDefinitions(artic);

ph = LMV.TRF.GetFeatureSet("phone")';
phIPA = MLing.ARPA2IPA(ph);

featTb = table(... 
    [repelem("articulatory", height(articTb)), repelem("phonemic", numel(ph))]', ...
    [articTb.fullName; phIPA], ...
    [articTb.shortName; ph], ...
    'VariableNames', ["Model", "Description", "Variable"]);

featTb.Description(11) = "relative pitch";
featTb.Description(12:22) = "rate of change in " + featTb.Description(1:11);
featTb.Description(23:end) = "binary impulse at the onset of " + featTb.Description(23:end);

featTb.Model(end+1) = "spectral";
featTb.Description(end) = "powers in Mel spectrogram in 80 bins from 0 to 8 kHz";
featTb.Variable(end) = "mel";

%% Save table

writetable(featTb, tbPath);
