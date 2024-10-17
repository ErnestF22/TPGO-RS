%function R=mean_rigid(allG,varargin)
%Compute the Karsher mean using gradient descent
%allG    4x4xN array of rigid body transformations
%Options
%weights Nx1 array of weights (default or [], all equal to 1)
%Nit     maximum number of iterations (default = Inf)
%init    initial estimation for the mean
function [G,A,B]=mean_rigid(allG,varargin)

Nit=Inf;
G=allG(:,:,1);
N=size(allG,3);
w=ones(N,1);

flagfixedNit=false;     %perform exactly Nit iterations (only if Nit<Inf)
flagnonmetric=false;    %use screw motions instead of the metric on SO(3) x R3

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'Nit'
            ivarargin=ivarargin+1;
            Nit=varargin{ivarargin};
        case 'init'
            ivarargin=ivarargin+1;
            G=varargin{ivarargin};
        case 'nonmetric'
            flagnonmetric=true;
        case 'weights'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
        case 'fixedNit'
            flagfixedNit=true;
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

prevAnorm=Inf;
it=0;
while it<Nit
    A=0;
    B=0;
    for iG=1:N
        if(w(iG)~=0)
            e=inv(G)*allG(:,:,iG);
            [a,b]=log_se3(e);
            if(~flagnonmetric)
                %b=e(1:3,4);
                b=allG(1:3,4,iG)-G(1:3,4);
            end
%             alla(:,iG)=a;
%             allb(:,iG)=b;
            A=A+w(iG)*a;
            B=B+w(iG)*b;
        end
    end
%     allb
    A=A/sum(w);
    B=B/sum(w);
    Anorm=norm([A;B]);
    if(Anorm>=prevAnorm & ~(flagfixedNit && Nit<Inf))
        break
    else
        if(~flagnonmetric)
            U=[rot(A) B; zeros(1,3) 1];
            %G=G*U;
            G=[G(1:3,1:3)*U(1:3,1:3) G(1:3,4)+U(1:3,4);zeros(1,3) 1];
        else
            U=exp_se3(A,B);
            G=G*U;
        end
        
        prevAnorm=Anorm;
    end
    it=it+1;
end
