%Given a stacked vector of variables, return only a named subset
function vSubstack=stackedVariables_substack(variables,nameCardinalityCell,v)
idxStacked=stackedVariables_positionsGet(variables,nameCardinalityCell);
vSubstack=v(idxStacked);
