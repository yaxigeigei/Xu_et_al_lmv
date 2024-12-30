classdef Artic
    %ARTIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        featNames = { ...
            'tt_x' 'tt_y' 'td_x' 'td_y' 'tb_x' 'tb_y' 'li_x' 'li_y' 'ul_x' 'ul_y' 'll_x' 'll_y' 'la' 'pro' 'ttcl' 'tbcl' 'v_x' 'v_y' ...
            'ja' 'ttcc' 'tbcc' 'tdcc'}'
    end
    
    methods(Static)
        function [defTb, dict] = GetFeatureDefinitions(arg)
            % Get the definitions of the estimated articulatory features
            % 
            %   [defTb, dict] = GetFeatureDefinitions()
            %   [defTb, dict] = GetFeatureDefinitions(featNames)
            %   [defTb, dict] = GetFeatureDefinitions(articTable)
            % 
            
            if ~exist('arg', 'var')
                featNames = NP.Artic.featNames;
            elseif istable(arg)
                featNames = arg.Properties.VariableNames;
            else
                featNames = arg;
            end
            featNames = string(featNames(:));
            
            dict.tt = 'tongue tip';
            dict.td = 'tongue dorsum';
            dict.tb = 'tongue body';
            dict.li = 'lower incisor';
            dict.ul = 'upper lip';
            dict.ll = 'lower tip';
            dict.la = 'lip aperture';
            dict.pro = 'lip protrusion';
            dict.ttcl = 'tongue tip constriction location';
            dict.tbcl = 'tongue body constriction location';
            dict.v = 'velum';
            dict.ja = 'jaw aperture';
            dict.ttcc = 'tongue tip constriction cosine';
            dict.tbcc = 'tongue body constriction cosine';
            dict.tdcc = 'tongue dorsum constriction cosine';
            
            defTb = table;
            defTb.shortName = featNames;
            for i = 1 : height(defTb)
                ss = strsplit(defTb.shortName(i), '_');
                if isfield(dict, ss{1})
                    ss(1) = string(dict.(ss{1}));
                else
                    ss(1) = defTb.shortName(i);
                    fprintf("%s is not a AKT variable.\n", defTb.shortName(i));
                end
                defTb.articulator(i) = ss(1);
                defTb.fullName(i) = strjoin(ss, ' ');
            end
        end
        
        function EnrichArticTable(se)
            % Derive AKT variables
            % 1) Additional track variables (TV) based on the following references
            %   Sivaraman et al. 2019. Unsupervised speaker adaptation for speaker independent acoustic to articulatory speech inversion.
            %   Parrot et al. 2019. Independent and automatic evaluation of acoustic-to-articulatory inversion models
            % 2) Time derivative from smoothed first-order quantities
            % 
            %   EnrichArticTable(se)
            % 
            
            tb = se.GetTable('artic');
            
            % Additional track variables (TV)
            tb.ja = cellfun(@NP.Artic.ComputeJA, tb.li_x, tb.li_y, tb.ul_x, tb.ul_y, 'Uni', false);
            tb.ttcc = cellfun(@NP.Artic.ComputeTXCC, tb.tt_x, tb.tt_y, 'Uni', false);
            tb.tbcc = cellfun(@NP.Artic.ComputeTXCC, tb.tb_x, tb.tb_y, 'Uni', false);
            tb.tdcc = cellfun(@NP.Artic.ComputeTXCC, tb.td_x, tb.td_y, 'Uni', false);
            
            
            % Remove previously computed derivatives
            isD = startsWith(tb.Properties.VariableNames, "d_");
            tb(:,isD) = [];
            
            % Compute time derivatives
            dtb = tb(:,2:end);
            dtb.Properties.VariableNames = "d_"+dtb.Properties.VariableNames;
            t = cat(1, tb.time{:});
            tSp = diff(t(1:2));
            order = 2;
            framelen = round(0.5/tSp/2)*2+1; % 0.05 sec; sgolayfilt requires even frame length
            
            for j = 1 : height(dtb)
                % Find ranges of AKT segments
                t = tb.time{j};
                m = diff(t) < tSp*2;
                bb = MMath.Logical2Bounds(m);
                bb(:,2) = bb(:,2)+1;
                
                % Piecewise smoothing
                for i = 1 : width(dtb)
                    v = double(dtb.(i){j});
                    for k = 1 : size(bb,1)
                        ind = bb(k,1):bb(k,2);
                        if numel(ind) >= framelen
                            v(ind) = sgolayfilt(v(ind), order, framelen);
                        end
                        v(ind) = gradient(v(ind)) ./ gradient(t(ind));
                    end
                    dtb.(i){j} = single(v);
                end
            end
            
            se.SetTable('artic', [tb dtb(:,2:end)]);
        end
        
        function ja = ComputeJA(li_x, li_y, ul_x, ul_y)
            % Compute jaw aperture from lower incisor and upper lip based on Sivaraman et al. 2019
            % "JA was defined as the Euclidean distance between the UL sensor and the LI sensor"
            % 
            %   ja = ComputeJA(li_x, li_y, ul_x, ul_y)
            % 
            ja = sqrt((li_x-ul_x).^2 + (li_y-ul_y).^2);
        end
        
        function txcc = ComputeTXCC(tx_x, tx_y)
            % Compute tongue constriction angle in cosine based on Parrot et al. 2019
            % e.g. "TTC: the cosine of the angle of the tongue tip off the horizontal axis"
            % 
            %   txcc = ComputeTXCC(tx_x, tx_y)
            % 
            txcc = tx_x ./ sqrt(tx_x.^2 + tx_y.^2);
        end
        
        function featLabels = GetLabels(featNames)
            % Convert feature names to standardized labels
            % 
            %   featLabels = GetLabels(featNames)
            % 
            
            featLabels = string(featNames);
            
            isFeat = startsWith(featNames, "d_");
            featLabels(isFeat) = erase(featNames(isFeat), "d_") + "'";
            
            isFeat = strcmp(featNames, 'drF0');
            featLabels(isFeat) = erase(featNames(isFeat), "d") + "'";
        end
        
    end
end

