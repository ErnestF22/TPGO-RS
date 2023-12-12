function [dw]=rotDyn_controlP_velocity(x,varargin)
flagCancelDynamics=false;
flagExtendedVector=false;
J=eye(3);
gainVelocityError=1;
flagVelocityInertiaWeight=false;
disturbance=zeros(12,size(x,2));

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'flagcanceldynamics'
            ivarargin=ivarargin+1;
            flagCancelDynamics=varargin{ivarargin};
        case 'flagvelocityinertiaweight'
            ivarargin=ivarargin+1;
            flagVelocityInertiaWeight=varargin{ivarargin};
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
            flagCancelDynamics=true;
        case 'gain'
            ivarargin=ivarargin+1;
            gainVelocityError=varargin{ivarargin};
        case 'gainvelocityerror'
            ivarargin=ivarargin+1;
            gainVelocityError=varargin{ivarargin};
        case 'flagextendedvector'
            ivarargin=ivarargin+1;
            flagExtendedVector=varargin{ivarargin};
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin}(1:12,:);
        otherwise
            % Skip data of unused entry
            ivarargin=ivarargin+1;
    end
    ivarargin=ivarargin+1;
end

[R,w]=rotDyn_stateUnpack(x);
if ~flagVelocityInertiaWeight
    uw=-w/norm(w);
else
    uw=-J*w/norm(w);
end
uw=gainVelocityError*uw;
dwDisturbance=rotDyn_inputUnpack(disturbance);
dw=uw-dwDisturbance;
if flagCancelDynamics
    dw=dw+rotDyn_gyroscopicTerm(w,'inertiaMatrix',J);
end

if flagExtendedVector
    dw=rotDyn_inputPack(dw,'state',R,w);
end

