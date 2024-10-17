%Given a cell array C of cell arrays, flattens it to a single cell array
%obtained by concatenating all the cell arrays
function CFlat=flattenTwoLevelsCell(C)
dims=cellfun(@(x) length(x),C);
CFlat=cell(1,sum(dims));
cumulativeDims=[0 cumsum(dims)];
for iC=1:length(C)
    CFlat(cumulativeDims(iC)+1:cumulativeDims(iC+1))=C{iC};
end

