% This function computed weighted unifrac distances between samples in
% relative_abundances
% relative_abundances is a matrix with relative abundances, where rows are samples, columns are taxa
% tree is, duh, the tree
% tree_info is the output of AAD_wUnifrac_get_tree_info. Have to input []
% for this function to get it if you don't provide it
% correction: 1 for using scaling factor as in Lozupone et al Applied and
% Env Micro 2007. 0 for not using it.
function out = AAD_wUnifrac(relative_abundances,tree,tree_info,correction)

% get tree info if not provided
if isempty(tree_info)
    tree_info = AAD_wUnifrac_get_tree_info(tree);
end

%
if logical(correction)
    d = nan(1,length(tree_info.leafnames));
    for i = 1:length(tree_info.leafnames)
        d(i) = sum(tree_info.distances(tree_info.upstreamnodes{i}));
    end
end

% calculate distances
nsamples = size(relative_abundances,1);

%
relative_abundances_all_branches=nan(nsamples,length(tree_info.leafnames)+length(tree_info.branch_leaves_idx));
for i = 1:nsamples
    sample = relative_abundances(i,:);
    totalbranchsample = nan(1,length(tree_info.branch_leaves_idx));
    for nbranch = 1:length(tree_info.branch_leaves_idx)
        totalbranchsample(nbranch) = sum(sample(tree_info.branch_leaves_idx{nbranch}));
    end
    relative_abundances_all_branches(i,:)=[sample,totalbranchsample];
end

out = nan(nsamples);
for i = 1:nsamples
    display(['working on row',num2str(i)])
    for j = i:nsamples
        difference = abs(relative_abundances_all_branches(i,:) - relative_abundances_all_branches(j,:));
        uncorrected_distance = full(sum(difference.*tree_info.distances'));
        if logical(correction)
            D = full(sum(d .* abs(sample1 - sample2)));
            out(i,j) = uncorrected_distance/D;
        else
            out(i,j) = uncorrected_distance;
        end
    end
end
for j=1:nsamples
    for i = j+1:nsamples
        out(i,j)=out(j,i);
    end
end



