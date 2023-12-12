%Builds the Jacobian of the constraints for a network
%function J=locSegJacobianConstraintsNetwork(E,x,R,constraintList,constraintParametersList)
%Inputs
%   E   [NEdges x 2] array with list of edges. Number of nodes is
%       determined as max(E(:)). Edges can be repeated. For constraints
%       that involve a single node, use E(iEdge,1)=E(iEdge,2)=iNode.
%   x   [2 x NNodes] or [3 x NNodes] array with locations of the nodes
%   constraintList  [1 x NEdges] cell-array of strings containing the kind of
%                   constraints for each edge in E
%   constraintParameters    [1 x NEdges] cell-array of cell-arrays containing additional
%                           parameters for the constraints. Each element is
%                           passed (if necessary) to the function computing
%                           the Jacobian for that constraint
%
% The Jacobian is computed at the location given by x and R
function J=locSegJacobianConstraintsNetwork(E,x,constraintList,constraintParametersList)
NEdges=size(E,1);
NNodes=max(E(:));
J=[];
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    xi=x(:,iNode);
    xj=x(:,jNode);

    constraintParameters=constraintParametersList{iEdge};
    Jr=locSegJacobianConstraintsBuild(constraintList{iEdge}, xi, xj, constraintParameters{:});
    JrExpanded=locSegJacobianConstraintsExpand(Jr,NNodes,iNode,jNode);
    J=[J; JrExpanded];
end
