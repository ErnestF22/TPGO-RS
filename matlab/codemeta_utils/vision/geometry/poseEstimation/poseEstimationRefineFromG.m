function GEst=poseEstimationRefineFromG(GInit,X,x,varargin)
optsMinimize={'maxIt',500};%,'showCost','displayIt'}; %
flagApproximatedHessian=true;
methodAbsolutePoses='pose';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optsminimize'
            ivarargin=ivarargin+1;
            optsMinimize=[optsMinimize{:} varargin{ivarargin}];
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if strcmpi(methodAbsolutePoses,'pose')
    GInit=invg(GInit);
end

f=@(G) cost(G,X,x,flagApproximatedHessian);
GEst=lie_minimizeGradNewton(rot3r3_funs(),f,GInit,optsMinimize{:});

if strcmpi(methodAbsolutePoses,'pose')
    GEst=invg(GEst);
end


function [c,gradc,Hessc]=cost(G,X,x,flagApproximatedHessian)
[c,gradc,Hessc]=poseEstimationCost_GFormat(G,X,x,'flagApproximatedHessian',flagApproximatedHessian);
c=sum(c);
gradc=sum(gradc,2);
Hessc=sum(Hessc,3);

