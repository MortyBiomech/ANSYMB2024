% ========================================================================
% FUNCTION: tfce2d - Threshold-Free Cluster Enhancement
% ------------------------------------------------------------------------
function TFCE = tfce2d(statMap, H, E, dh)
    % Vectorised 2-D TFCE implementation (Smith & Nichols, 2009)
    
    if nargin<4, dh = 0.1; end
    if nargin<3, E = 0.5; end
    if nargin<2, H = 2; end
    
    statMap(statMap<0) = 0;             % TFCE defined for positive stats
    hVals = 0:dh:max(statMap(:));
    TFCE = zeros(size(statMap));
    
    % 8-connected neighbourhood in 2D
    conn = conndef(2,'maximal');
    
    for h = hVals
        mask = statMap >= h;            % Threshold at height h
        CC = bwconncomp(mask, conn);    % Identify clusters
        
        for clIdx = 1:CC.NumObjects
            voxIdx = CC.PixelIdxList{clIdx};
            extent = numel(voxIdx);
            TFCE(voxIdx) = TFCE(voxIdx) + (h^H) * (extent^E) * dh;
        end
    end
end