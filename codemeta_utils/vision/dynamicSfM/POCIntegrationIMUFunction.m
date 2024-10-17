
function [R,T]=POCIntegrationIMUFunction(t,w,alphaIMU,R0,T0,v0,varargin)
flagRight=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'left'
            flagRight=false;
        case 'right'
            flagRight=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

g=gravityVector();

R=rotDyn_integrateVelocity(t,w,R0,'flagRight',flagRight);
a=multiprodMatVec(invR(R),alphaIMU)-repmat(g,1,size(alphaIMU,2));
T=realDyn_integrateAcceleration(t,a,T0,v0);
