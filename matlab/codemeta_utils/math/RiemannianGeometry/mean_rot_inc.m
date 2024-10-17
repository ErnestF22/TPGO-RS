%function R=mean_rot(allR,varargin)
%Compute the Karsher mean using gradient descent
%allR    3x3xN array of rotations
%Options
%weights Nx1 array of weights (default or [], all equal to 1)
%Nit     maximum number of iterations (default = Inf)
%init    initial estimation for the mean
function [R,A]=mean_rot_inc(allR, varargin)

Nit=Inf;
R=allR(:,:,1);
N=size(allR,3);
w=ones(N,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'Nit'
            ivarargin=ivarargin+1;
            Nit=varargin{ivarargin};
        case 'init'
            ivarargin=ivarargin+1;
            R=varargin{ivarargin};
        case 'weights'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

prevAnorm=Inf;
it=0;
S=eye(3);
while it<Nit
    A=0;
    for iR=1:N
        if(w(iR)~=0)
            A=A+w(iR)*logrot(S'*R'*allR(:,:,iR));
        end
    end
    A=A/sum(w);
    Anorm=norm(A);
    if(Anorm>=prevAnorm)
        break
    else
        S=S*rot(A);
        prevAnorm=Anorm;
    end
    it=it+1;
end
R=R*S;
