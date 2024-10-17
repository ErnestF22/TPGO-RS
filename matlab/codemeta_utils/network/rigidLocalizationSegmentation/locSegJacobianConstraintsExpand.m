%Splits a Jacobian matrix and puts the two parts in a row-block
%function JrowBlock=locSegJacobianConstraintsExpand(Jr,N,iNode,jNode)
function JrowBlock=locSegJacobianConstraintsExpand(Jr,N,iNode,jNode)
K=size(Jr,1);
D=size(Jr,2)/2;

JrowBlock=zeros(K,N*D);
JrowBlock(:,D*(iNode-1)+(1:D))=Jr(:,1:D);
if jNode~=iNode
    JrowBlock(:,D*(jNode-1)+(1:D))=Jr(:,(D+1):end);
end
