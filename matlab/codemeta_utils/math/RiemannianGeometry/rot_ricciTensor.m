%Evaluates the Ricci tensor Ric(X,Y)
%function r=rot_ricciTensor(R,X,Y)
%Inputs
%   R       the rotation at which the tangent vectors are based
%   X,Y     two tangent vectors
%If X and Y are not supplied, returns the matrix representation of the
%tensor in the basis given by rot_tangentBasis(R)
function r=rot_ricciTensor(R,X,Y)
T=rot_tangentBasis(R);
dim=size(T,3);

if exist('X','var')
    r=ricciTensorEval(R,X,Y,T);
else
    r=zeros(dim);
    for iDim=1:dim
        Ei=T(:,:,iDim);
        for jDim=1:dim
            Ej=T(:,:,jDim);
            r(iDim,jDim)=ricciTensorEval(R,Ei,Ej,T);
        end
    end
end

function r=ricciTensorEval(R,X,Y,T)
dim=size(T,3);
r=0;
for kDim=1:dim
    Ek=T(:,:,kDim);
    r=r+rot_curvature(R,X,Ek,Y,Ek);
end
