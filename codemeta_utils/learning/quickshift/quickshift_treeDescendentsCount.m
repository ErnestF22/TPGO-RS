%Count number of descendents for each datapoint
%function cnt=quickshift_treeDescendentsCount(descendents)
function cnt=quickshift_treeDescendentsCount(descendents)
cnt=cellfun(@length,descendents,'UniformOutput',true);
