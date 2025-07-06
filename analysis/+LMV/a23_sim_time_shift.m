%% Plot RNN simulation results

anaDir = LMV.Data.GetAnalysisDir("rnn");

%% Load simulated response shift data

respTb = readtable(fullfile(anaDir, "data", "Response_time_trajectories_1s.csv"));

%% Plot response shifts

cycles = unique(respTb.cycle);
nCycles = numel(cycles);

f = MPlot.Figure(1326); clf
tl = tiledlayout(nCycles, 1);
tl.Padding = "compact";

for i = 1 : nCycles
    ax = nexttile;
    
    cyc = cycles(i);
    mask = respTb.cycle == cyc;
    t = respTb.time_rel_s(mask);
    y = respTb.Rmu(mask);
    
    % Plot Rmu vs time
    plot(ax, t, y, 'k'); hold on;
    
    % Draw solid vertical bar at t_peak_s (alpha 0.3)
    t_peak = unique(respTb.t_peak_s(mask));
    yl = ylim(ax);
    plot(ax, [t_peak t_peak], yl, '-', 'Color', [0 0 0 .3]);
    
    % Set axis limits
    ax.XLim = [-.2 .2];
    ax.YLim = [0 500];
    
    % Add trial number as text in the top left
    xText = ax.XLim(1) + 0.01;
    yText = ax.YLim(2) - 0.05 * diff(ax.YLim);
    text(ax, xText, yText, sprintf('Trial %d', cyc), 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    ax.YLabel.String = num2str(i+1);
    ax.YLabel.Rotation = 0;
    ax.YTickLabel = [];
    if i == nCycles
        ax.XLabel.String = 'Time from input onset (s)';
    else
        ax.XLabel.String = '';
        ax.XTickLabel = [];
        ax.XAxis.Visible = 'off';
    end
    MPlot.Axes(ax);
end

MPlot.Paperize(f, .33, .75);
exportgraphics(f, fullfile(anaDir, "response_shifts.png"));
exportgraphics(f, fullfile(anaDir, "response_shifts.pdf"), "ContentType", "vector");
