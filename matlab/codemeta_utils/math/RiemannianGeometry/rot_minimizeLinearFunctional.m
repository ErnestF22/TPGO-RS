%Minimize linear functional on O(n)
%function rot_minimizeLinearFunctional(A,b,R0,varargin)
%Minimizes a functional of the form norm(A*R(:)-b,2)
function [R,output]=rot_minimizeLinearFunctional(A,b,R0,varargin)
d=size(R0,1);
EHat=reshape(rot_tangentBasis(eye(d)),d^2,[]);
[R,output]=lie_minimizeGradNewton(rot_funs(),@(R) cost(R,A,b,EHat),...
    R0,varargin{:});

function [c,g,H]=cost(R,A,b,EHat)
f=A*R(:)-b;
c=(f'*f)/2;
if nargout>1
    g=-EHat'*kron(eye(3),R')*A'*f;
    if nargout>2
        H=hessian(R,A,f,EHat);
    end
end

function H=hessian(R,A,f,EHat)
n=size(R,1);
d=size(EHat,2);
H=zeros(d,d);
for k=1:d
    H(:,k)=EHat'*kron(eye(3),reshape(EHat(:,k),n,n)'*R')*A'*f;
end
H=H+EHat'*(kron(eye(3),R')*(A'*A)*kron(eye(3),R))*EHat;
