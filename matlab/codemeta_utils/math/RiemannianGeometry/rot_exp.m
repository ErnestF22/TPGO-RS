%function R=rot_exp(R1,A)
%Compute exponential map of the tangent vector A at R1 
%Note: A is assumed to be of the form R1*skew-symmetric matrix

%%AUTORIGHTS%%

function R=rot_exp(R1,A,varargin)

%these flags control tricks to improve the accuracy
flagRotMult=true;           %use numerically stable but slow (3x slower) rotation multiplication
flagEnforceSkew=true;       %enforce tangent vector to be skew symmetric
flagAverageSVDBasis=false;   %average U and V in SVD decomposition of the tangent vector
flagOrthSVDBasis=true;     %make sure V in the SVD decomposition is orthogonal

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagrotmult'
            ivarargin=ivarargin+1;
            flagRotMult=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N2=size(A,3);
if(N2>1)
    %recursive call if A contains multiple tangent vectors
    if size(R1,3)==1
        R=zeros([size(R1) N2]);
        for iN2=1:N2
            R(:,:,iN2)=rot_exp(R1,A(:,:,iN2));
        end
    else
        R=zeros(size(R1));
        for iN2=1:N2
            R(:,:,iN2)=rot_exp(R1(:,:,iN2),A(:,:,iN2));
        end
    end       
else
    N=size(A,1);
    R=zeros(size(A));

    if ~isempty(R1)
        A=R1'*A;
    end
    
    if flagEnforceSkew
        A=(A-A')/2;
    end
    
    [U,S,V]=svd(A);
    
    if flagAverageSVDBasis && norm(A,'fro')>1e-15
        for is=1:2:N-1
            s=sign(U(:,is)'*V(:,is+1));
            V(:,is)=(V(:,is)-s*U(:,is+1))/2;
            V(:,is+1)=(V(:,is+1)+s*U(:,is))/2;
        end
        if mod(N,2)==1
            s=sign(U(:,end)'*V(:,end));
            V(:,end)=(V(:,end)+s*U(:,end))/2;
        end
    end
    
    if flagOrthSVDBasis
        [UV,SV,VV]=svd(V);
        V=UV*VV';
    end
    
    for is=1:2:N-1
        s=S(is,is)*sign(U(:,is)'*V(:,is+1));
        ctheta=cos(s);
        stheta=sin(s);
        R(is:is+1,is:is+1)=[ctheta -stheta; stheta ctheta];
    end
    
    if mod(N,2)==1
        R(end,end)=1;
    end
    
    if flagRotMult
        Rt=R';              %workaround for a known bug in Matlab R2009b
        Vt=V';
        RV=(Rt\V'+R/V)/2;
        R=(Vt\RV+V/RV')/2;
    else
        R=V*R*V';
    end
    
    if ~isempty(R1)
        if flagRotMult
            R1t=R1';
            R=(R1t\R+R1/R')/2;
        else
            R=R1*R;
        end
    end
end
