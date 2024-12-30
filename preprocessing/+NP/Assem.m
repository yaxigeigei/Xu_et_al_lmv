classdef Assem
    
    methods(Static)
        function [A, meta] = ComputeAssembly(R)
            %
            
            % z-score spike rates
            Z = MMath.Normalize(R, 'zscoresoft', NP.Param.normAddMax);
            
            % PCA
            sPCA = struct;
            sPCA.subNames = cellstr("PC" + (1:size(Z,2)));
            [sPCA.B, PCproj, sPCA.eigVal, ~, sPCA.varExplained] = pca(Z, 'Centered', false);
            sPCA.varExplained = sPCA.varExplained(:);
            
            % Determine the number of components using Marcenko-Pastur threshold
            [nT, nU] = size(Z);
            th = (1 + sqrt(nU/nT))^2;
            nC = sum(sPCA.eigVal > th);
            
            % Compute rICA
            Zp = PCproj(:, 1:nC);
            Mdl = rica(Zp, nC);
            
            % Compute neuronal weights for assemblies
            V = sPCA.B(:,1:nC) * Mdl.TransformWeights;
            
            % Normalize weight vectors
            % columns of V are already unit vectors
            
            % % Unify signs such that
            % [~, I] = max(abs(W));
            % ind = sub2ind(size(W), I, 1:nC);
            
            % Determine assembly membership by a threshold
            thV = std(V)*2 + mean(abs(V));
            isV = abs(V) > thV;
            
            % Compute assembly activation
            A = zeros(nT, nC);
            parfor k = 1 : nC
                
                P = V(:,k) * V(:,k)';
                P = P - diag(diag(P));
                
                for i = 1 : nT
                    A(i,k) = Z(i,:) * P * Z(i,:)';
                end
            end
            
            % 
            meta = struct;
            meta.PCA = sPCA;
            meta.ICA.Mdl = Mdl;
            meta.ICA.V = V;
            meta.ICA.thV = thV;
            meta.ICA.isV = isV;
        end
        
        
    end
end

