function [n,v,Av,An,b]=homFlowParametersEstimate4ptTranslationRefine(n,v,x,dx,w,varargin)
maxIt=20;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NX=size(x,2);
preAv=[];
preAn=[];
b=[];
for iX=1:NX
    xi=x(1,iX);
    yi=x(2,iX);
    A2i=[eye(2) -x(:,iX)]';
    A1i=[x(:,iX)' 1];
    bi=dx(:,iX)-[-w(2)+yi*w(3)+xi*yi*w(1)-xi^2*w(2);
        w(1)-xi*w(3)+yi^2*w(1)-xi*yi*w(2)];
    
    preAn=[preAn;kron(A1i',A2i')];
    preAv=[preAv;kron(vec(A2i'),A1i)];
    b=[b;bi];
end

for it=1:maxIt
    AnVec=reshape(preAn*v,2,3,[]);
    An=reshape(permute(AnVec,[1 3 2]),[],3);
    n=An\b;
    
    AvVec=reshape(preAv*n,2,3,[]);
    Av=reshape(permute(AvVec,[1 3 2]),[],3);
    v=Av\b;
end
