%% Plot synaptic variables and efficacy from CSV data

anaDir = LMV.Data.GetAnalysisDir("rnn");

%% Load synaptic variable data

synTb = readtable(fullfile(anaDir, "data", "synaptic_variables_data.csv"));

%% Plot synaptic variables

% Extract variables
T = synTb.time_s;
T = T - T(1);
Rmu = synTb.Rmu;
u = synTb.u;
x = synTb.x;
A = synTb.A;
uxA = synTb.uxA;

% Create figure with 3 vertical panels
f = MPlot.Figure(1326); clf
tl = tiledlayout("vertical");
tl.Padding = 'compact';

% --- Top panel: Pre-synaptic input ---
ax1 = nexttile;
plot(ax1, T, Rmu, 'Color', [0 0.2 0.6]);
ax1.YLabel.String = 'R_\mu (Hz)';
ax1.YLabel.FontSize = 18;
ax1.YLabel.Color = [0 0.4 0.8];
ax1.YLabel.FontAngle = 'italic';
ax1.XTickLabel = [];
ax1.Title.String = 'Pre-synaptic input';
ax1.Title.FontSize = 16;
ax1.Title.FontWeight = 'normal';
ax1.Box = 'off';
ax1.XLim = [min(T) max(T)];
ax1.YLim = [0 max(Rmu)*1.1];

% --- Middle panel: Synaptic variables ---
ax2 = nexttile;

yyaxis(ax2, 'left');
plot(ax2, T, u, 'Color', [0.2 0.8 0.2]); hold on; % green
plot(ax2, T, x, 'Color', [1 0 0], 'LineStyle', '-');                % red, solid
ax2.YColor = [0 0 0]; % axis ticks and numbers in black
ax2.YLabel.String = '\color{green}u, \color{red}x';
ax2.YLabel.FontSize = 16;
ax2.YLabel.FontAngle = 'italic';
ax2.YLim = [0 1];

yyaxis(ax2, 'right');
plot(ax2, T, A, 'Color', [0.2 0.4 1]);            % blue
ax2.YColor = [0.2 0.4 1]; % axis ticks and numbers in black
ax2.YLabel.String = 'A';
ax2.YLabel.Color = [0.2 0.4 1];
ax2.YLabel.FontSize = 16;
ax2.YLabel.FontAngle = 'italic';
ax2.YLim = [0 40];

ax2.Title.String = 'Synaptic variables';
ax2.Title.FontSize = 16;
ax2.Title.FontWeight = 'normal';
ax2.XTickLabel = [];
ax2.Box = 'off';
ax2.XLim = [min(T) max(T)];

% --- Bottom panel: Synaptic efficacy ---
ax3 = nexttile;
plot(ax3, T, uxA, 'k');
ax3.YLabel.String = 'uxA';
ax3.YLabel.FontSize = 16;
ax3.YLabel.FontAngle = 'italic';
ax3.XLabel.String = 'Time (s)';
ax3.XLabel.FontSize = 18;
ax3.Title.String = 'Synaptic efficacy';
ax3.Title.FontSize = 16;
ax3.Title.FontWeight = 'normal';
ax3.Box = 'off';
ax3.XLim = [min(T) max(T)];
ax3.YLim = [0 max(uxA)*1.1];

% Adjust font sizes for all axes
axs = [ax1 ax2 ax3];
for ax = axs
    MPlot.Axes(ax);
end

MPlot.Paperize(f, .45, .45);
exportgraphics(f, fullfile(anaDir, "synaptic_vars.png"));
exportgraphics(f, fullfile(anaDir, "synaptic_vars.pdf"), "ContentType", "vector");
