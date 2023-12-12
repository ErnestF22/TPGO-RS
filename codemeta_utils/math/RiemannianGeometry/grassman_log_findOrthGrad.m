function [YY2,errors]=grassman_log_findOrthGrad(YY1,YY2,p,varargin)
Nit=100;
t=1.1;
tol=1e-15;

flagOptimalLineSearch=false;
flagCollectErrors=false;
flagLinearInit=false;
flagRotMult=true;              %optional argument for rot_exp and rot_log

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optimallinesearch'
            flagOptimalLineSearch=true;
        case 'linearinit'
            flagLinearInit=true;
        case 'stepsize'
            ivarargin=ivarargin+1;
            t=varargin{ivarargin};
        case 'nit'
            ivarargin=ivarargin+1;
            Nit=varargin{ivarargin};
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


if(nargout>1)
    flagCollectErrors=true;
end
    
n=size(YY1,1);
cost=@(YY2) norm(projVert(rot_log([],YY1'*YY2,'flagRotMult',flagRotMult),p))^2;

if(flagLinearInit)
    %Find a good initialization as solution to a Procrustes problem
    [U,~,V]=svd(YY1(:,1:p)'*YY1(:,1:p));
    R1=U*V';
    [U,~,V]=svd(YY2(:,p+1:end)'*YY1(:,p+1:end));
    R2=U*V';
    R=blkdiag(R1,R2);
    %if(det(R)<0)
        %[U,S,V]=svd(U*diag([ones(p-1,1);-1])*S*V');
        %R=[U*V'];
    %end
    if(det(R)>0)
        YY2=YY2*R;
    end
end
    

if(det(YY2)<0)
    %this should not happen
    error('YY2 has negative determinant, check the code!')
end

In=eye(n);

for it=1:Nit
    %Project gradient on the vertical space
    Hg=projVert(rot_log([],YY1'*YY2,'flagRotMult',flagRotMult),p);
    if(norm(imag(Hg),'fro')>0)
        error('Hg has imaginary part!')
    end
    
    if(flagOptimalLineSearch)
        [t,costval]=fminbnd(@(t) cost(YY2*rot_exp(In,-t*Hg,'flagRotMult',flagRotMult)),0,2);
        YY2=YY2*rot_exp(In,-t*Hg);
    else
        YY2=YY2*rot_exp(In,-t*Hg,'flagRotMult',flagRotMult);
        costval=cost(YY2);
    end
    
    if(flagCollectErrors)
        errors.allcost(it)=costval;
        errors.allnorm(it)=norm(t*Hg)^2;
    end
    %disp(costval)

    if(costval<tol^2)
        break;
    end
    
end

function H=projVert(H,p)
H(p+1:end,1:p)=0;
H(1:p,p+1:end)=0;

