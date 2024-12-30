classdef Param
    
    properties(Constant)
        % Stimulus
        stimDict = struct( ...
            "mdlc2_si2244", "nobody likes snakes", ...
            "msjs1_si1899", "you took me by surprise", ...
            "fcaj0_si1479", "have you got enough blankets", ... % x2
            "mfxv0_si1635", "we've got plenty of time to think about that", ... % x2
            "fsjk1_si2285", "the girl nodded understandingly", ...
            "fcag0_si1641", "were you in love with that girl", ...
            "mafm0_si2199", "will you tell me why", ...
            "mbbr0_si2315", "junior what on earth's the matter with you", ... % x2
            "mewm0_si1978", "something pulled my leg", ...
            "fdxw0_si2141", "it was nobody's fault", ...
            "fltm0_si2330", "i'm going to search this house", ...
            "mrgg0_si1199", "they close sometime after eight", ...
            "mbom0_si1644", "they even pay me six dollars a month", ... % x2
            "madd0_si1295", "he may get a tax refund", ...
            ... % (The following stim were removed since NP41)
            "mrew1_si2130", "does this bother you", ...
            "fsrh0_si1719", "we've done our part", ...
            "mjmm0_si625",  "a tiny handful never did make the concert", ... % x2
            "fgmd0_si2107", "it sounded silly why go on", ...
            "mdmt0_si1832", "you're right now buddy", ...
            "mdpk0_si552",  "we think differently", ...
            "mkcl0_si1721", "come home right away", ...
            "mclk0_si1660", "the car door crashed shut", ...
            "mmdm0_si1941", "we'll talk over at your office", ...
            ... % (The following stim were removed since NP38)
            "mtjg0_si2157", "who would take over", ...
            "mmgk0_si1952", "we're going someplace", ...
            "mjac0_si2148", "look somewhere else", ...
            "mrlj1_si2332", "i saw your horse outside", ...
            "mrgs0_si1986", "what's your name, anyway", ...
            "mtxs0_si1690", "hey, come back, he shouted" ...
            );
        
        stimIdList = string(fieldnames(LMV.Param.stimDict));
        stimIdList14 = LMV.Param.stimIdList(1:14);
        stimIdList4 = LMV.Param.stimIdList([3 4 8 13]);
        stimIdList12 = LMV.Param.stimIdList([3 4 8 13 1 2 5:7 9:11]);
        
        stimTextList = string(struct2cell(LMV.Param.stimDict));
        stimTextList14 = LMV.Param.stimTextList(1:14);
        stimTextList4 = LMV.Param.stimTextList([3 4 8 13]);
        stimTextList12 = LMV.Param.stimTextList([3 4 8 13 1 2 5:7 9:11]);
        
        % Responsiveness
        respBaselineDur = 0.3; % in sec
        respWinOffset = 0.2; % in sec
        respWinMinDur = 0.2; % in sec
        
        % 
        regions = ["mPrCG", "vPrCG", "IFG", "STG"];
        phases = ["atten"];
    end
    
    methods(Static)
        function [uId, uIdUni] = GetSelectedClusId(recId)
            % Return handpicked cluster IDs for a given recording
            % 
            %   [uId, uIdUni] = GetSelectedClusId(recId)
            % 
            % Input
            %   recId           Recording ID such as 'NP41_B1'.
            % Outputs
            %   uId             Short original cluster IDs from Kilosort/Phy/MTracer.
            %   uIdUni          Long and globally unique cluster IDs with recording and probe digits added.
            % 
            switch recId
                case 'NP35_B1', uId = [106 179 72 48]; % v2
                case 'NP35_B2', uId = [898 349 906 909 792 795 814 815 850 110 705 487 838 931]; % v2
                case 'NP38_B5', uId = [266 21 220 296 66 25 202 189 2 278 277]; % v2
                case 'NP38_B6', uId = [11 105 103 96 89 157 45 177 82 9 75]; % v2
                case 'NP41_B1', uId = [245 233 224 205 202 168 155 138 134 434 441 446 48 455 16 409 367 418 188 171 463 450 432 454 64]; % v1.5
                case 'NP41_B2', uId = []; % TBD. This recording only has 32 trials
                case 'NP42_B1', uId = []; % v2
                case 'NP43_B1', uId = [210 161 99 570 403 578 579 598 613 133 619 27 428 413 582 583 584 154 131 624 544 319 227 552 560 612 549 307]; % v2
                case 'NP44_B2', uId = [173 137 115 479 287 284 492 248 501 507 528 532 462 473 484 523]; % v2
                case 'NP44_B3', uId = [515 295 236 574 575 585 587 594 267 601 604 619 602 620]; % v2
                case 'NP45_B1', uId = [254 269 276 322 91 286 65 298 205 165 323 314 48]; % v2
                case 'NP45_B2', uId = [307 608 287 283 613 617 52 583 597]; % v2
                case 'NP46_B1', uId = [145 128 129 240 98 30 13 247 3 250 248 0 219 242 258 236 263 269 271 246 251 283]; % v2
                case 'NP47_B3', uId = [4 133 102 99 62 139 19 145 12]; % v1.5
                case 'NP50_B3', uId = [74 265 452 446 454 456 458 459 464]; % v2
                case 'NP51_B1', uId = []; % no clear response patterns
                case 'NP52_B1', uId = [94 45 40 101 16 6 103 3 87 51 56 105 34 89 49]; % v2?
                case 'NP52_B2', uId = [71 150 152 153 154 0]; % v2?
                case 'NP53_B1', uId = [35 32 17 13 11 195 189]; % v2_mx
                case 'NP54_B1', uId = [153 78 54 34 28 161 135 129 167 168 99 1]; % v2; auto [44 34 27 16 108 89 87 79 3 70 43 0]
                case 'NP55_B1', uId = [106 170 110 85 53 167 143 131 13 207 49 75 70 59 114 99 86 189 206]; % auto
                case 'NP55_B2', uId = []; % auto
                case 'NP56_B1', uId = [[888 897 1049 160 948 426 1067 390 402 855], ... % v2
                        1e4 + [615 618 633 640 644 351 658 683 163 706 714 726 737 510 295 747 748 826 789 450 607 75]]; % v2
                case 'NP60_B1', uId = [393 722 537 884 912 282 938 575 564 899 920 845]; % v2
                case 'NP69_B1', uId = [[490 316 242 147 443 404 279 249 519 511], ... % v2_mx
                        1e4 + [739 98 746 67 73 750 91 774 26 189 737 773 197]]; ... % v2_mx
                case 'NP69_B2', uId = [[152 32 342 228 214 374 163 164 292 159], ... % v2_mx
                        1e4 + [206 12 22 27 26 240 175 149 114 56 53 447 436 58 400]]; ... % v2_mx
                case 'NP74_B1', uId = [[61 160 348 170 44 205 50 454 470 337 7 6], 1e4 + [460 238 10 336 412]]; % auto
                case 'NP78_B1', uId = [44 17 69 57]; % v2_qg
                case 'NP81_B1', uId = []; % v2_mx
                case 'NP94_B1', uId = [10199 10083 10028 10394]; % auto
                % case 'NP106_B1', uId = [328 613 290 179 174, 1e4+[383 296 643 14 735 183], 2e4+[0]]; % lmv, auto
                case 'NP106_B1', uId = [445 599 291 230 220 177 145 52 31 677 646 507 182, ...
                        1e4+[189 77 36 643 14 634 626 183 603 736 308], ...
                        2e4+[0 365 333 214 195]]; % lmv-br, auto
                case 'NP122_B1', uId = [193 189 176 271 131 126 128 127 266 273]; % auto
                otherwise, uId = [];
            end
            uIdUni = NP.Unit.GetBaseClusId(recId) + uId;
        end
        
        function uid = GetExampleClusId(anaName)
            % Return handpicked cluster IDs for a given analysis
            % 
            %   uid = GetExampleClusId(anaName)
            % 
            % Input
            %   anaName     Name of an analysis.
            % Output
            %   uid         Long, globally unique cluster IDs including the recording and probe digits.
            % 
            switch anaName
                case 'selectivity'
%                     uid = [ ...
%                         460100240, ... % 0-1mm
%                         460100098 440200462 450100286 460100013 460100005, ... % 1-2mm
%                         410100249 410100260 440200479 460100003 410100245 460100248, ... % 2-3mm
%                         440300575 410100211 440300585 440300575 410100367 450200307, ... % 3-4mm
%                         410100463 410100441 410100170 450100323 440300587, ... % 4-5mm
%                         410100441 410100454 440300601 410100434 470300139, ... % 5-6mm
%                         ];
                    uid = [ ...
                        410100245 450200307 460100248 460100030 460100005 450100286, ...  % phasic units. 
                        410100155 410100441 440300575 410100450 450100323 410100463, ...  % non-phasic units. 
                        ];
                case 'sen4'
                    
                otherwise
                    uid = [];
            end
        end
        
        function cc = GetRegionColors(regions)
            % Get n-by-3 color codes for brain regions using "lines" colormap.
            % 
            %   cc = GetRegionColors(regions)
            % 
            %   mPrCG: yellow
            %   vPrCG: orange
            %   IFG: purple
            %   STG: blue
            % 
            regions = cellstr(regions);
            cmap = lines();
            for i = numel(regions) : -1 : 1
                switch regions{i}
                    case 'mPrCG'
                        cc(i,:) = cmap(3,:); % yellow
                    case 'vPrCG'
                        cc(i,:) = cmap(2,:); % orange
                    case 'IFG'
                        cc(i,:) = cmap(4,:); % purple
                    case 'STG'
                        cc(i,:) = cmap(1,:); % blue
                    otherwise
                        cc(i,:) = [0 0 0];
                end
            end
        end
        
        function cc = GetTaskPhaseColors(n)
            % Get n-by-3 color codes for task phases using Dark2 colormap.
            % 
            %   cc = GetTaskPhaseColors(numPhases)
            %   cc = GetTaskPhaseColors(phaseNames)
            % 
            %   1. stim: blueish green
            %   2. delay: yellow
            %   3. init: brown
            %   4. prod: magenta
            %   5. iti: violet
            %   6. (unassigned): gray
            % 
            
            % Make colormap
            cmap = brewermap([], 'Dark2');
            
            if isnumeric(n)
                % Get a given number of colors in order
                cc = repmat(cmap, ceil(n/size(cmap,1)), 1);
                cc = cc(1:n,:);
            else
                % Get colors by phase names
                n = string(n);
                for i = numel(n) : -1 : 1
                    switch n(i)
                        case 'stim'
                            cc(i,:) = cmap(1,:);
                        case 'delay'
                            cc(i,:) = cmap(6,:);
                        case 'init'
                            cc(i,:) = cmap(7,:);
                        case 'prod'
                            cc(i,:) = cmap(4,:);
                        case {'iti', 'atten'}
                            cc(i,:) = cmap(3,:);
                        otherwise
                            cc(i,:) = cmap(8,:);
                    end
                end
            end
        end
        
        function cc = GetSentenceColors(n)
            % Get n-by-3 color codes for different sentences using Set3 colormap.
            % 
            %   cc = GetSentenceColors(numSen)
            %   cc = GetSentenceColors(stimIds)
            % 
            
            % Make colormap
            stimIdList = LMV.Param.stimIdList14;
            cmap = brewermap(12, 'Paired');
            cmap = [cmap; mean(cmap(1:2,:)); mean(cmap(3:4,:)); 0.7 0.7 0.7];
            
            if isnumeric(n)
                % Get a given number of colors in order
                cc = repmat(cmap, ceil(n/size(cmap,1)), 1);
                cc = cc(1:n,:);
            else
                % Get colors by stimId
                stimIds = string(n);
                cInd = zeros(size(n));
                for i = 1 : numel(cInd)
                    k = find(stimIds(i)==stimIdList, 1);
                    if isempty(k)
                        k = numel(stimIdList)+1;
                    end
                    cInd(i) = k;
                end
                cc = cmap(cInd,:);
            end
        end
        
        function cc = GetLinkerColors(types)
            % Get n-by-3 color codes for linker types using Accent colormap.
            % 
            %   cc = GetLinkerColors(types)
            % 
            %   1. mirror: green
            %   2. bridge: orange
            %   3. feedback: purple
            %   6. (unassigned): gray [.9 .9 .9]
            % 
            types = cellstr(types);
            cmap = brewermap([], 'Accent');
            for i = numel(types) : -1 : 1
                switch types{i}
                    case 'mirror'
                        cc(i,:) = cmap(1,:); % green
                    case 'bridge'
                        cc(i,:) = cmap(3,:); % orange
                    case 'feedback'
                        cc(i,:) = cmap(2,:); % purple
                    otherwise
                        cc(i,:) = [0 0 0] + 0.5;
                end
            end
        end
        
        function cc = GetModelColors(mdlNames)
            % Get n-by-3 color codes for linker types using Set2 colormap.
            % 
            %   cc = GetModelColors(mdlNames)
            % 
            %   Articulatory models:    pink
            %   Acoustic models:        green
            %   Phoneme models:         brown
            %   Other models:           gray
            % 
            mdlNames = cellstr(mdlNames);
            cmap = brewermap([], 'Set2');
            for i = numel(mdlNames) : -1 : 1
                switch mdlNames{i}
                    case {'artic', 'artic3'}
                        cc(i,:) = cmap(2,:);
                    case 'strf'
                        cc(i,:) = cmap(5,:);
                    case 'phone'
                        cc(i,:) = cmap(7,:);
                    otherwise
                        cc(i,:) = [0 0 0] + 0.5;
                end
            end
        end
        
    end
end
