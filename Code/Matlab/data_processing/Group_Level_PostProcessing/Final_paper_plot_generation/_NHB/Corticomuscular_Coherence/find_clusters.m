function [pixels_per_cluster, masses] = find_clusters(binary_map, value_map, conn)
%   binary_map : [nF x nT] logical, true where pixel is above threshold
%   value_map  : [nF x nT] real, the values to sum within each cluster
%   conn       : 4 or 8 (pixel connectivity)
%
%   Returns:
%     pixels_per_cluster : cell array, linear indices of each cluster
%     masses             : vector of cluster masses (sum of value_map)

CC = bwconncomp(binary_map, conn);
pixels_per_cluster = CC.PixelIdxList;
nC = numel(pixels_per_cluster);

masses = zeros(nC, 1);
for c = 1:nC
    masses(c) = sum(value_map(pixels_per_cluster{c}));
end
end