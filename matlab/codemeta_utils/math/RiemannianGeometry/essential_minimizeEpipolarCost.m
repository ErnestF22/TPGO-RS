function QEst=essential_minimizeEpipolarCost(Q0,x1,x2,varargin)
optsMinimize={'maxIt',500,'stopcriterion','normGrad',...
    'stopThreshold',1e-15};%,'showLogCost','displayIt'};
flagApproximatedHessian=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optsminimize'
            ivarargin=ivarargin+1;
            optsMinimize=[optsMinimize{:} varargin{ivarargin}];
        case 'flagapproximatedhessian'
            ivarargin=ivarargin+1;
            flagApproximatedHessian=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagApproximatedHessian
    f=@(Q) cost(Q,x1,x2);
else
    f=@(Q) costApproximatedHessian(Q,x1,x2);
end    
QEst=lie_minimizeGradNewton(essential_funs(),f,Q0,optsMinimize{:});

function [c,gradc,Hessc]=cost(Q,x1,x2)
[c,gradc,Hessc]=essential_evaluateEpipolarCost(Q,x1,x2);
c=sum(c);
gradc=sum(gradc,2);
Hessc=sum(Hessc,3);

function [c,gradc,Hessc]=costApproximatedHessian(Q,x1,x2)
[c,gradc,Hessc]=essential_evaluateEpipolarCost(Q,x1,x2,'approximatedHessian');
c=sum(c);
gradc=sum(gradc,2);
Hessc=sum(Hessc,3);
