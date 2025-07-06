%% Plot RNN simulation results

anaDir = LMV.Data.GetAnalysisDir("rnn");

%% Load simulated trial data

% Load activity traces
tsTb = readtable(fullfile(anaDir, "data", "activity_traces.csv"));

% Load on/off times
evtTb = readtable(fullfile(anaDir, "data", "on_off_times.csv"));

%% Plot timeseries

% Define cluster order: persistent first, then activity-silent
clusOrder = [1 3 2 4]; % MATLAB is 1-indexed, so Cluster0=1, Cluster1=2, etc.
clusNames = {'Persistent activity neuron 1', ...
                 'Persistent activity neuron 3', ...
                 'Activity-silent neuron 2', ...
                 'Activity-silent neuron 4'};
clusColors = [
    0.8 0.2 0.2;   % red-ish for cluster 0
    0.2 0.4 0.8;   % blue-ish for cluster 1
    1.0 0.6 0.4;   % orange-ish for cluster 2
    0.2 0.7 1.0    % cyan-ish for cluster 3
    ];
ylims = [0 100; 0 100; 0 200; 0 200];

f = MPlot.Figure(1324); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

for i = 1 : numel(clusOrder)
    ax = nexttile;
    clusIdx = clusOrder(i);
    clusCol = sprintf('Cluster%d', clusIdx-1); % -1 for zero-indexed in CSV
    
    % Plot activity trace
    x = tsTb.time_s;
    y = tsTb.(clusCol);
    plot(x, y, 'Color', clusColors(clusIdx,:), 'LineWidth', 1.5);
    hold on;

    yl = ylim; % Get current y-limits for line height

    for k = 1:height(evtTb)
        cidx = evtTb.cluster(k) + 1; % cluster is 0-indexed in table, 1-indexed in MATLAB
        this_color = clusColors(cidx, :);
        % t_load
        plot([evtTb.t_load_s(k) evtTb.t_load_s(k)], [0 yl(2)], '--', 'Color', this_color, 'LineWidth', 1.5);
        % t_unload
        plot([evtTb.t_unload_s(k) evtTb.t_unload_s(k)], [0 yl(2)], '--', 'Color', this_color, 'LineWidth', 1.5);
    end
    
    % Set y-limits: lower = 0, upper = auto
    ylim(ylims(i,:));
    
    % Label subplot
    ylabel('Ru (Hz)');
    if i == 4
        xlabel('Time (s)');
    else
        ax.XLabel = [];
        ax.XTickLabel = [];
    end
    ax.XLim = [1 6.4];
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 0.8, 0.6);
exportgraphics(f, fullfile(anaDir, "activity_traces.png"));
exportgraphics(f, fullfile(anaDir, "activity_traces.pdf"), "ContentType", "vector");
