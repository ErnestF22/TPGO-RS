%Compute function of residuals and derivatives of the rotation-only trifocal constraint
%function [f,gf,Dgf]=trifocalRotationCost(R1,R2,n0,n1,n2,funs)
%Compute the cost funs.f(n0'*hat(R1*n1)*R2*n2), its gradient and its
%differential of the gradient w.r.t. R1,R2.
%Inputs
%   n0,n1,n2    [3xN] arrays of N normals in \real{3}
%   R1,R2       [3x3] rotations
%   funs        Structure with three fields: f, df and ddf, which are
%               handles to a function, its first and second derivatives.
%               By default, funs.f=@(t) 0.5*t^2;
%Outputs
%   f           [1xN] array where f(i)=funs.f(n0(:,i)'*hat(R1*n1(:,i))*R2*n2(:,i))
%   gf          [6xN] array where fe(:,i) is the gradient of f(i)
%   Dgf         [6x6xN] array where Dfe(:,:,i) is the differential of fe(:,i)
function [f,gf,Dgf]=trifocalRotationCost(R1,R2,n0,n1,n2,funs)
if ~exist('funs','var')
    funs.f=@(t) t.^2/2;
    funs.df=@(t) t;
    funs.ddf=@(t) ones(size(t));
end

Nn=size(n0,2);
flagComputeJacobian=false;
flagComputeSecondJacobian=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end


if ~flagComputeJacobian
    e=trifocalRotationResiduals(R1,R2,n0,n1,n2);
else
    if ~flagComputeSecondJacobian
        [e,ge]=trifocalRotationResiduals(R1,R2,n0,n1,n2);
    else
        [e,ge,Dge]=trifocalRotationResiduals(R1,R2,n0,n1,n2);
        Dgf=zeros(size(Dge));
    end
    df=zeros(size(e));
    gf=zeros(size(ge));
end
f=zeros(size(e));

for in=1:Nn
    f(in)=funs.f(e(in));
    if flagComputeJacobian
        df(in)=funs.df(e(in));
        gf(:,in)=df(in)*ge(:,in);
        if flagComputeSecondJacobian
            Dgf(:,:,in)=funs.ddf(e(in))*ge(:,in)*ge(:,in)'+df(in)*Dge(:,:,in);
        end
    end
end
