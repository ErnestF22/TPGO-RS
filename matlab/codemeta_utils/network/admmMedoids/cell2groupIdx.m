%Given a cell array of arrays, return a cell array of the form {[1 ... 1],
%[2 .... 2],...}, where the length of each individual array is equal to the
%number of columns in the corresponding element of the initial cell array.
function xi0Group=cell2groupIdx(xi0)
nbGroups=length(xi0);
sz=cellfun(@(x) size(x,2), xi0);
xi0Group=cell(nbGroups);
for iGroup=1:nbGroups
    xi0Group{iGroup}=iGroup*ones(1,sz(iGroup));
end
end
