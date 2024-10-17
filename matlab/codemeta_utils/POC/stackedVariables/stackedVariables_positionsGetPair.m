%Give indeces in the outer product of stacked variables
%function [idxStacked,idxStackedTranspose]=stackedVariables_positionsGetPair(variables,nameCardinalityCellRow,nameCardinalityCellCol)
%Let M=v*v' be the matrix with outer product of stacked variables.
%idxStacked returns the matrix of linear indeces identifying the entries in
%M corresponding to the rows and columns of the selected variables
function [idxStacked,idxStackedTranspose]=stackedVariables_positionsGetPair(variables,nameCardinalityCellRow,nameCardinalityCellColumn)
d=stackedVariables_stackSize(variables);

idxStackedRow=stackedVariables_positionsGet(variables,nameCardinalityCellRow);
idxStackedCol=stackedVariables_positionsGet(variables,nameCardinalityCellColumn);

idxStacked=idxStackedRow+d*(idxStackedCol-1)';
idxStackedTranspose=idxStackedCol+d*(idxStackedRow-1)';
