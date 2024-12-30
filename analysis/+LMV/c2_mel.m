%% 

srcTb = LMV.Data.FindSource('NP46_B1');
se = NP.SE.LoadSession(srcTb.path{1});

%% 

NP.Audio.AddMelTable(se);

%% 

[tt, tv] = se.GetTable("taskTime", "taskValue");
tr = find(tv.stimText=="we've got plenty of time to think about that", 1);
tWin = [tt.stimOn(tr) tt.stimOff(tr)];

f = MPlot.Figure(123); clf
ax = nexttile;
NP.TaskBaseClass.PlotMelSpectrogram(ax, se, tr, tWin, 'Channel', "speaker1", 'Wave', 0);
% ax.XLabel.String = [];
% ax.XTickLabel = [];
MPlot.Paperize(f, 1.5, 0.3);

