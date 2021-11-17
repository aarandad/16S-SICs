% This function gets the relevant information about a tree to calculate
% weighted Unifrac
% The only input is the tree
% out is a structure with fields:
% - upstreamnodes: for each node, what are the nodes that take it up to the
% root of the tree. 
% - branch_leaves_idx: leaf index for all leaves in a given branch
% - branch_leaves_names: leaf name for all leaves in a given branch
% - leafnames: leaf names in tree. Can also get using get(tree,'leafnames')
% - distances: distances for each branch. Can also get using
% get(tree,'distances')
function out = AAD_wUnifrac_get_tree_info(tree)

leafnames = get(tree,'LeafNames');
pointers = get(tree,'Pointers');
branchnames = get(tree,'BranchNames');
nodenames = get(tree,'NodeNames');
distances = get(tree,'distances');
upstreamnodes = cell(1,length(leafnames));
for i = 1:length(leafnames)
    idxtemp = find(strcmp(nodenames,leafnames{i}));
    idxbranch = find(pointers(:,1)==idxtemp|pointers(:,2)==idxtemp);
    upstreamnodes_temp = idxtemp;
    while ~isempty(idxbranch)
        idxtemp =find(strcmp(branchnames{idxbranch},nodenames));
        upstreamnodes_temp = [upstreamnodes_temp,idxtemp];
        idxbranch = find(pointers(:,1)==idxtemp|pointers(:,2)==idxtemp);
    end
    upstreamnodes{i} = upstreamnodes_temp;
end

branch_leaves_names = cell(1,length(branchnames));
branch_leaves_idx = cell(1,length(branchnames));
for i = 1:length(branchnames)
    tr = subtree(tree,getbyname(tree,branchnames{i}));
    branch_leaves_names{i} = get(tr,'LeafNames');
    for j = 1:length(branch_leaves_names{i})
        branch_leaves_idx{i}(j) = find(strcmp(branch_leaves_names{i}{j},leafnames));
    end
end


out.upstreamnodes = upstreamnodes';
out.branch_leaves_idx = branch_leaves_idx';
out.branch_leaves_names = branch_leaves_names';
out.leafnames = leafnames;
out.distances = distances;






