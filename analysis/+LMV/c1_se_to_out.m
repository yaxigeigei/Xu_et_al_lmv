%% Save out data for Matt

srcTb = NP.Data.FindSource('lmv', 'Region', ["mPrCG", "STG"]);

%% Convert original se

outDir = "C:\chang_lab\project_np\preproc\out_structs\original";
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for k = 4 %1 : height(srcTb)
    % Check output
    recId = srcTb.recId{k};
    outPath = fullfile(outDir, recId+"_out.mat");
    if exist(outPath, 'file')
        fprintf("\nThe out struct for %s already exists.\n", recId);
        continue
    end
    
    % Load se
    se = NP.SE.LoadSession(srcTb.path{k});
    
    % Conversion
    out = LMV.SE.ToOutStruct(se);
    
    % Save out
    save(outPath, '-v7.3', 'out');
end

%% Convert time morphed se

outDir = "C:\chang_lab\project_np\preproc\out_structs\morphed";
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for k = 1 : height(srcTb)
    % Check output
    recId = srcTb.recId{k};
    outPath = fullfile(outDir, recId+"_out_m1.mat");
    if exist(outPath, 'file')
        fprintf("\nThe out struct for %s already exists.\n", recId);
        continue
    end
    
    % Load se
    sePath = fullfile(NP.Data.GetAnalysisRoot, "data", "se_m1", recId+"_se_m1.mat");
    se = NP.SE.LoadSession(sePath);
    
    % Conversion
    out = LMV.SE.ToOutStruct2(se);
    
    % Save out
    save(outPath, '-v7.3', 'out');
end

return
%% 

r = arrayfun(@(x) x.spikes_sua(9,:)', out(2:end), 'Uni', false)';
prodOn = arrayfun(@(x) find(x.prodOn,1), out(2:end))';
stimOn = arrayfun(@(x) find(x.stimOn,1), out(2:end))';

dur = prodOn - stimOn;
[~, I] = sort(dur);
r = r(I);

f = MPlot.Figure(123); clf
r0 = (1:numel(r))*25;
MPlot.PlotTraceLadder([], r, r0, 'ColorArray', 'k'); hold on
plot(stimOn(I), r0, 'b', 'LineWidth', 1.5)
plot(prodOn(I), r0, 'r', 'LineWidth', 1.5)
MPlot.Axes(gca);
axis tight





