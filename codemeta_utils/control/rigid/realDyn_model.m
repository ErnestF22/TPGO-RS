%Implement the second order model for a point mass
%function dx=realDyn_model(x,varargin)
%The model is given by m*ddx=F+Delta
function dx=realDyn_model(x,varargin)
F=zeros(3,1);
m=1;
disturbance=zeros(6,size(x,2));

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'mass'
            ivarargin=ivarargin+1;
            m=varargin{ivarargin};
        case 'force'
            ivarargin=ivarargin+1;
            F=varargin{ivarargin};
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NPoses=size(x,2);
if size(F,2)==1 && NPoses>1
    F=repmat(F,1,NPoses);
end
[~,v]=realDyn_stateUnpack(x);
[vDisturbance,FDisturbance]=realDyn_stateUnpack(disturbance);
dx=[v+vDisturbance;(F+FDisturbance)/m];
