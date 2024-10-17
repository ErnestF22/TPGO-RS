%Apply matrix(ces) to vector(s)
%function Mv=multiprodMatVec(M,v)
%Computes Mv(:,i)=M(:,:,i)*v(:,i). If size(M,3)==1 or size(v,2)==1, these
%are automatically expanded to match the dimension of the other argument.
function Mv=multiprodMatVec(M,v)
Mv=multiprod(M,v,[1 2],1);
