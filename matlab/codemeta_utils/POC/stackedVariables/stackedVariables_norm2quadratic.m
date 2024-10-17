%Transform squared L2 norm of residuals into quadratic expression
%function [Q,p,c]=stackedVariables_norm2quadratic(variables,names,A,b,type)
%Transform norm(A*vNames-b)^2 (type=='point',default) or
%norm(vNames'*A-b)^2 (type='model') into v'*Q*v+2*p'*v+c, where vNames is a
%vector containing only the variables listed in names, and v is a vector
%with all variables
function [Q,p,c]=stackedVariables_norm2quadratic(variables,nameCardinalitiesCell,A,b,type)
if ~exist('type','var') || isempty(type)
    type='point';
end
[QReduced,pReduced,c]=norm2quadratic(A,b,type);
dStack=stackedVariables_stackSize(variables);
Q=zeros(dStack);
p=zeros(dStack,1);
idxStack=stackedVariables_positionsGet(variables,nameCardinalitiesCell);
Q(idxStack,idxStack)=QReduced;
p(idxStack)=pReduced;


