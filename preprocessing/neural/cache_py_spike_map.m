
function cache_py_spike_map(ksFolder, rez)
    % Save spike map and motion correction data for plotting in Python
    % 
    %   cache_py_spike_map(ksDir)
    %   cache_py_spike_map(ksDir, rez)
    % 
    % Inputs
    %   ksDir           The path of Kilosort output folder that contains the "rez.mat" file.
    %   rez             The rez struct. If not provided, it will be loaded from "rez.mat".

    % Load data
    if nargin < 2
        load(fullfile(ksFolder, "rez.mat"), "rez");
    end

    % Prepare spike map data
    tb = table();
    if isfield(rez, 'st0')
        st = rez.st0;
        tb.time = (st(:,1) + rez.ops.tstart) / rez.ops.fs;
        tb.y = st(:,2);
        tb.amp = st(:,3);
    else
        fprintf("Raw drift map data (rez.st0) is not available. Plot sorted spike map (rez.st3).\n");
        st = rez.st3;
        tb.time = st(:,1) / ops.fs;
        iTemp = st(:,2);
        iChan = rez.iNeighPC(16, iTemp);
        tb.y = rez.ycoords(iChan);
        tb.amp = round(st(:,3));
    end
    
    % Save spike map as a readable table
    tb = sortrows(tb, {'time', 'y'});
    outputFile = fullfile(ksFolder, 'spike_map.csv');
    writetable(tb, outputFile);

    % Save motion correction traces
    if ~isfield(rez, 'dshift')
        fprintf("Motion correction data (rez.dshift) is not available.\n");
    else
        writeNPY(rez.dshift, fullfile(ksFolder, 'motion.npy'));
    end
    
    fprintf("Done! Results save in %s\n", ksFolder);
end
