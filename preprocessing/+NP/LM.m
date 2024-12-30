classdef LM
    methods(Static)
        function coords = ScaleCoords(coords, subjectID, expId)
            % Convert from pixels to millimeters
            switch subjectID
                case 'NP00'
                    k = 1;
                otherwise
                    k = 1;
                    warning('This recording does not have a defined scaling factor');
            end
            coords = coords / k;
        end
        
        function s = GetInd()
            % Get the MediaPipe FaceMesh indices for the landmarks of interest
            
            % Indices go from subject's left to right
            % Lip arc # goes from inside to outside
            s.upLipArc4 = [292 410 271 270 268   1  38  40  41 186  62]';
%             s.upLipArc3 = [307 409 305 304 303  12  73  74  75 185  77]';
%             s.upLipArc2 = [293 408 273 272 269  13  39  42  43 184  63]';
            s.upLipArc1 = [309 416 311 312 313  14  83  82  81 192  79]';
            s.loLipArc1 = [309 325 319 403 318  15  88 179  89  96  79]';
%             s.loLipArc2 = [293 326 320 404 317  16  87 180  90  97  79]';
%             s.loLipArc3 = [307 308 321 405 316  17  86 181  91  78  77]';
            s.loLipArc4 = [292 376 322 406 315  18  85 182  92 147  62]';
            
            % Indices go from subject's left to right
            s.chinArc = [366 380 379 401 378 153 149 177 150 151 137]';
            
            % Index at nose tip
            s.noseTip = 5;
            
%             % Indices go from lateral to medial
%             s.upLeftEye = [264 467 389 388 387 386 385 399 363]';
%             s.loLeftEye = [264 250 291 374 375 381 382 383 363]';
%             s.upRightEye = [34 247 162 161 160 159 158 174 134]';
%             s.loRightEye = [34   8 164 145 146 154 155 156 134]';
            
            % Indices go in ascending order
            s.leftIris = (475 : 478)';
            s.rightIris = (470 : 473)';
        end
        
        function sXY = FindXY(coords)
            % Find XY coordinates for the landmarks of interest
            coords = permute(coords, [3 1 2]); % change dims to time-by-points-by-xyz
            sInd = NP.LM.GetInd();
            fnList = fieldnames(sInd);
            for i = 1 : numel(fnList)
                fn = fnList{i};
                sXY.([fn 'X']) = coords(:, sInd.(fn), 1);
                sXY.([fn 'Y']) = coords(:, sInd.(fn), 2);
            end
        end
        
        function [mouthHeight, mouthWidth] = GetMouthSize(sXY)
            % Get the height and width of the mouth opening
            %
            %   [mouthHeight, mouthWidth] = GetMouthSize(sXY)
            % 
            
            % Find relevant coordinates
            upMidXY = [sXY.upLipArc1X(:,6) sXY.upLipArc1Y(:,6)];
            loMidXY = [sXY.loLipArc1X(:,6) sXY.loLipArc1Y(:,6)];
            leftXY = [sXY.upLipArc1X(:,1) sXY.upLipArc1Y(:,1)];
            rightXY = [sXY.upLipArc1X(:,end) sXY.upLipArc1Y(:,end)];
            
            % Vertical distance between the two inner midpoints
            vVect = upMidXY - loMidXY;
            for i = size(vVect,1) : -1 : 1
                mouthHeight(i,1) = norm(vVect(i,:));
            end
            
            % Horizontal distance between the two inner corner points
            hVect = rightXY - leftXY;
            for i = size(hVect,1) : -1 : 1
                mouthWidth(i,1) = norm(hVect(i,:));
            end
        end
        
        function d = GetChin2NoseDist(sXY)
            % Get the distance from nose tip to mid chin
            %
            %   d = GetChin2NoseDist(sXY)
            % 
            chinMidXY = [sXY.chinArcX(:,6) sXY.chinArcY(:,6)];
            noseTipXY = [sXY.noseTipX sXY.noseTipY];
            v = chinMidXY - noseTipXY;
            for i = size(v,1) : -1 : 1
                d(i,1) = norm(v(i,:));
            end
        end
        
        function [leftX, leftY, rightX, rightY] = GetPupilPositions(sXY)
            % Get the X and Y coordinates of left and right pupils
            %
            %   [leftX, leftY, rightX, rightY] = GetPupilPositions(sXY)
            % 
            leftX = mean(sXY.leftIrisX, 2);
            leftY = mean(sXY.leftIrisY, 2);
            rightX = mean(sXY.rightIrisX, 2);
            rightY = mean(sXY.rightIrisY, 2);
        end
        
        function AddManualPupil(se, folderPath)
            % Add manually tracked pupil location to 
            
            lm = se.GetTable('landmark');
            placeholder = cellfun(@(x) NaN(size(x)), lm.time, 'Uni', false);
            lm.pupilX = placeholder;
            lm.pupilY = placeholder;
            lm.pupilV = placeholder;
            
            lmSearch = MBrowse.Dir2Table(fullfile(folderPath, "* landmarks.mat"));
            for i = 1 : height(lmSearch)
                % Read files
                load(fullfile(lmSearch.folder{i}, lmSearch.name{i}), "landmarksTable");
                pp = landmarksTable.pupil;
                for j = 1 : numel(pp)
                    if isempty(pp{j})
                        pp{j} = [NaN NaN];
                    end
                end
                pp = cat(1, pp{:});
                
                % Get trial index
                nameParts = strsplit(lmSearch.name{i}, '_');
                k = str2double(nameParts{1}(6:end)); % after the letters of 'trial'
                
                % Add to epoch
                lm.pupilX{k} = pp(:,1);
                lm.pupilY{k} = pp(:,2);
                
                % Derive velocity
                ppv = [0 0; diff(pp)];
                fs = se.userData.landmarkMeta.frame_rate;
                lm.pupilV{k} = sqrt(ppv(:,1).^2 + ppv(:,2).^2) * fs; % in pixels/sec
            end
            
            se.SetTable('landmark', lm);
        end
        
        % Plots
        function PlotMeasurements(ax, sXY)
            % 
            
            [mouthH, mouthW] = NP.LM.GetMouthSize(sXY);
            chin2nose = NP.LM.GetChin2NoseDist(sXY);
            t = (0 : numel(mouthH)-1)' * 1/29.97;
            
            cla(ax)
            plot(ax, t, [mouthH mouthW chin2nose]);
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Distance (px)';
            MPlot.Axes(ax);
        end
        
        function PlotLandmark(ax, se, trialIdx, t, tfObj)
            % 
            
            % Get landmark data at the requested timepoint
            tb = se.GetTable('landmark');
            s = table2struct(tb(trialIdx,:));
            [~, idx] = min(abs(s.time - t));
            edgeLen = double(se.userData.landmarkMeta.img_size(1));
            s = structfun(@(x) x(idx,:) * edgeLen, s, 'Uni', false);
            
            % Contour of lips and chin
            arcNames = {'upLipArc1', 'upLipArc4', 'loLipArc1', 'loLipArc4', 'chinArc'};
            for i = numel(arcNames) : -1 : 1
                xx{i} = s.(arcNames{i}+"X")';
                yy{i} = s.(arcNames{i}+"Y")';
            end
            xx = cat(2, xx{:});
            yy = cat(2, yy{:});
            
%             % Contour of iris
%             irNames = {'leftIris', 'rightIris'};
%             xxIris = [s.leftIrisX([1:end 1]) NaN s.rightIrisX([1:end 1])]';
%             yyIris = [s.leftIrisY([1:end 1]) NaN s.rightIrisY([1:end 1])]';
%             xx = [xx xxIris];
%             yy = [yy yyIris];
            
            % Height and width of mouth opening
            xMouH = [s.upLipArc1X(6) s.loLipArc1X(6)]';
            yMouH = [s.upLipArc1Y(6) s.loLipArc1Y(6)]';
            xMouW = [s.upLipArc1X(1) s.loLipArc1X(end)]';
            yMouW = [s.upLipArc1Y(1) s.loLipArc1Y(end)]';
            
            % Nose to chin dist
            xChin = [s.noseTipX s.chinArcX(6)]';
            yChin = [s.noseTipY s.chinArcY(6)]';
            
            % Pupils
            xPu = [s.leftPupilX s.rightPupilX]';
            yPu = [s.leftPupilY s.rightPupilY]';
            
            % Transform
            if exist('tfObj', 'var') && ~isempty(tfObj)
                [xx, yy] = transformPointsForward(tfObj, xx, yy);
                [xMouH, yMouH] = transformPointsForward(tfObj, xMouH, yMouH);
                [xMouW, yMouW] = transformPointsForward(tfObj, xMouW, yMouW);
                [xChin, yChin] = transformPointsForward(tfObj, xChin, yChin);
                [xPu, yPu] = transformPointsForward(tfObj, xPu, yPu);
            end
            
            % Overwrite with manual pupil position
            xPu = [s.pupilX s.pupilX]' / edgeLen;
            yPu = [s.pupilY s.pupilY]' / edgeLen;
            if ~isfield(ax.UserData, 'puTfObj')
                [~, ax.UserData.puTfObj] = Img23.Transform(zeros(1280), 'Translate', [-175 250], 'Crop', [300 300]);
            end
            puTfObj = ax.UserData.puTfObj;
            [xPu, yPu] = transformPointsInverse(puTfObj, xPu, yPu);
            if exist('tfObj', 'var') && ~isempty(tfObj)
                [xPu, yPu] = transformPointsForward(tfObj, xPu, yPu);
            end
            
            % Update plot
            if isfield(ax.UserData, 'vert') && ~isempty(ax.UserData.vert) && ishandle(ax.UserData.vert)
                for i = 1 : size(xx,2)
                    ax.UserData.arcs(i).XData = xx(:,i);
                    ax.UserData.arcs(i).YData = yy(:,i);
                end
                ax.UserData.vert.XData = xMouH;
                ax.UserData.vert.YData = yMouH;
                ax.UserData.hori.XData = xMouW;
                ax.UserData.hori.YData = yMouW;
                ax.UserData.chin.XData = xChin;
                ax.UserData.chin.YData = yChin;
                ax.UserData.pu.XData = xPu;
                ax.UserData.pu.YData = yPu;
            else
                ax.UserData.arcs = plot(ax, xx, yy, 'w', 'LineWidth', 1);
                cc = lines(10);
                ax.UserData.vert = plot(ax, xMouH, yMouH, 'o', 'Color', cc(1,:), 'LineWidth', 2);
                ax.UserData.hori = plot(ax, xMouW, yMouW, 'o', 'Color', cc(2,:), 'LineWidth', 2);
                ax.UserData.chin = plot(ax, xChin, yChin, 'o', 'Color', cc(3,:), 'LineWidth', 2);
                ax.UserData.pu = plot(ax, xPu, yPu, 'o', 'Color', cc(5,:), 'LineWidth', 2);
            end
        end
        
        function PlotLandmarkFromCoords(ax, coords)
            % 
            
            k = round(ax.UserData.time * 29.97);
            if k < 1 || k > size(coords,3)
                return
            end
            c = coords(:, 1:2, k);
            sXY = NP.LM.FindXY(c);
            
            % Contour of lips and chin
            xx = [sXY.upLipArc1X' sXY.upLipArc4X' sXY.loLipArc1X' sXY.loLipArc4X' sXY.chinArcX'];
            yy = [sXY.upLipArc1Y' sXY.upLipArc4Y' sXY.loLipArc1Y' sXY.loLipArc4Y' sXY.chinArcY'];
            
            % Height and width of mouth opening
            xMouH = [sXY.upLipArc1X(6) sXY.loLipArc1X(6)]';
            yMouH = [sXY.upLipArc1Y(6) sXY.loLipArc1Y(6)]';
            xMouW = [sXY.upLipArc1X(1) sXY.loLipArc1X(end)]';
            yMouW = [sXY.upLipArc1Y(1) sXY.loLipArc1Y(end)]';
            
            % Nose to chin dist
            xChin = [sXY.noseTipX sXY.chinArcX(6)]';
            yChin = [sXY.noseTipY sXY.chinArcY(6)]';
            
            % Update plot
            if isfield(ax.UserData, 'vert') && ~isempty(ax.UserData.vert) && ishandle(ax.UserData.vert)
                for i = 1 : size(xx,2)
                    ax.UserData.arcs(i).XData = xx(:,i);
                    ax.UserData.arcs(i).YData = yy(:,i);
                end
                ax.UserData.vert.XData = xMouH;
                ax.UserData.vert.YData = yMouH;
                ax.UserData.hori.XData = xMouW;
                ax.UserData.hori.YData = yMouW;
                ax.UserData.chin.XData = xChin;
                ax.UserData.chin.YData = yChin;
            else
                ax.UserData.arcs = plot(ax, xx, yy, 'w', 'LineWidth', 1);
                
                cc = lines(3);
                ax.UserData.vert = plot(ax, xMouH, yMouH, 'Color', cc(1,:), 'LineWidth', 2);
                ax.UserData.hori = plot(ax, xMouW, yMouW, 'Color', cc(2,:), 'LineWidth', 2);
                ax.UserData.chin = plot(ax, xChin, yChin, '--', 'Color', cc(3,:), 'LineWidth', 2);
                
                ax.XLim = [min(c(:,1)), max(c(:,1))] + [-1 1]*50;
%                 ax.YLim = [min(c(:,2)), max(c(:,2))] + [-1 1]*50;
                ax.YLim = [600 1200];
            end
        end
        
    end
    
end


