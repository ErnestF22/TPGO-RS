%function x=linsysorth(A,b,x0)
%Solves Ax=b subject to x'*x=1 iteratively
function x=linsysorth(A,b,varargin)

x0=A\b;
epsilon=0.1;

flagshowplots=false;    %show plots of lambda and x
display='off';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'epsilon'
            ivarargin=ivarargin+1;
            epsilon=varargin{ivarargin};
        case 'init'
            ivarargin=ivarargin+1;
            x0=varargin{ivarargin};
        case 'showplots'
            flagshowplots=true;
        case 'display'
            ivarargin=ivarargin+1;
            display=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Optional argument not valid!')
    end
    ivarargin=ivarargin+1;
end

AtA=A'*A;
Atb=A'*b;
x=x0;
%lambda=(x'*A'*b-x'*AtA*x)/(x'*x);
lambda=0;

I=eye(size(AtA));

x = fmincon(@costx,x,[],[],[],[],[],[],@nonlinconst,optimset('MaxFunEvals',1e6,'LargeScale','off','Display',display,'TolX',1e-6));

% lambda=fminunc(@cost,lambda,optimset('GradObj','on','Display','off'));
% x=(AtA+lambda*I)\Atb;

% for(it=1:1000)
%     lambda=lambda+epsilon*(x'*x-1);
%     x=(AtA+lambda*I)\Atb;
%     
%     if(flagshowplots)
%         alllambda(it)=lambda;
%         allx(:,it)=x;
%     end
% end
% 
% if(flagshowplots)
%     figure(1)
%     plot(alllambda)
% 
%     figure(2)
%     plot(allx')
% end


    function [f,g]=cost(lambda)
    x=(AtA+lambda*I)\Atb;
    f=-0.5*(sum((A*x-b).^2)+lambda*(x'*x-1));
    if(nargout>1)
        g=-(x'*x-1);
    end
    end

    function f=costx(x)
        f=sum((A*x-b).^2);
    end

    function [c, ceq]=nonlinconst(x)
        c=0;
        ceq=x'*x-1;
    end
end
