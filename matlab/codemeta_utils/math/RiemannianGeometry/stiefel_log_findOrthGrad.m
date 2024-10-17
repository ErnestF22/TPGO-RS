function [YY2,errors]=stiefel_log_findOrthGrad(YY1,YY2,p,varargin)
Nit=100;
t=1.1;
tol=1e-12;

flagOptimalLineSearch=false;
flagCollectErrors=false;
flagLinearInit=true;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optimallinesearch'
            flagOptimalLineSearch=true;
        case 'stepsize'
            ivarargin=ivarargin+1;
            t=varargin{ivarargin};
        case 'nit'
            ivarargin=ivarargin+1;
            Nit=varargin{ivarargin};
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        case 'nolinearinit'
            flagLinearInit=false;
        case 'linearinit'
            flagLinearInit=true;
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
cost=@(YY2) norm(projVert(rot_log([],YY1'*YY2),p))^2;

if(flagLinearInit)
    %Find a good initialization as solution to a Procrustes problem
    [U,~,V]=svd(YY2(:,p+1:end)'*YY1(:,p+1:end));
    R=U*V';
    %if(det(R)<0)
        %[U,S,V]=svd(U*diag([ones(p-1,1);-1])*S*V');
        %R=[U*V'];
    %end
    if(det(R)>0)
        YY2=YY2*blkdiag(eye(p),R);
    end
end

if(det(YY2)<0)
    %this should not happen
    error('YY2 has negative determinant, check the code!')
end

Ip=eye(p);
Inp=eye(n-p);

for it=1:Nit
    %Project gradient on the vertical space
    %Hg=projVert(rot_log([],YY1'*YY2),p);
    Hg=rot_log([],YY1'*YY2);
    if(norm(imag(Hg),'fro')>0)
        error('Hg has imaginary part!')
    end
    Hg=Hg(p+1:end,p+1:end);
    
    if(flagOptimalLineSearch)
        [t,costval]=fminbnd(@(t) cost(YY2*blkdiag(Ip,rot_exp(Inp,-t*Hg))),0,2);
        YY2=YY2*blkdiag(Ip,rot_exp(Inp,-t*Hg));
    else
        
        YY2(p+1:end,p+1:end)=YY2(p+1:end,p+1:end)*rot_exp(Inp,-t*Hg); %optimized version of YY2=YY2*blkdiag(Ip,rot_exp(Inp,-t*Hg));
        costval=cost(YY2);
    end
    
    if(flagCollectErrors)
        errors.allcost(it)=costval;
        errors.allnorm(it)=norm(Hg)^2;
    end

    if(sqrt(costval)<tol)
        break;
    end
    
end

function H=projVert(H,p)
H(:,1:p)=0;
H(1:p,:)=0;
