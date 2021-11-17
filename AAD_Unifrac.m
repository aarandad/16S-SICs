% This function computed weighted unifrac distances between samples in
% relative_abundances
% relative_abundances is a matrix with relative abundances, where rows are samples, columns are taxa
% tree is, duh, the tree
% tree_info is the output of AAD_wUnifrac_get_tree_info. Have to input []
% for this function to get it if you don't provide it
function out = AAD_Unifrac(relative_abundances,tree,tree_info)

% get tree info if not provided
if isempty(tree_info)
    tree_info = AAD_wUnifrac_get_tree_info(tree);
end

nsamples = size(relative_abundances,1);

% get nodes 
presentnodes = cell(nsamples,1);
for i = 1:nsamples
    sample = relative_abundances(i,:);
    idx = find(sample>0);
    presentnodes{i} = unique(cell2mat(tree_info.upstreamnodes(idx)'));
    
end
out = nan(nsamples);
for i = 1:nsamples
    display(['working on sample ',num2str(i)])
    for j = i:nsamples
        nodes_idx1 = presentnodes{i};
        nodes_idx2 = presentnodes{j};
      
        idx_nodes_unique = setxor(nodes_idx1,nodes_idx2);
        idx_nodes_observed = unique([nodes_idx1,nodes_idx2]);
        
        out(i,j) = sum(tree_info.distances(idx_nodes_unique))./sum(tree_info.distances(idx_nodes_observed));
        
        
    end
end


for j=1:nsamples
    for i = j+1:nsamples
        out(i,j)=out(j,i);
    end
end