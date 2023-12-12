function [dw]=rotDyn_controlPD_contraction(x,varargin)
% PD controller inside contraction region
flagCancelDynamics=false;
flagExtendedVector=false;
RReference=eye(3);
wReference = zeros(3,1);
J=eye(3);
gainRotationError=1;
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
        case 'rreference'
            ivarargin=ivarargin+1;
            RReference=varargin{ivarargin};
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
            flagCancelDynamics=true;
        case 'gain'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
            gainVelocityError=varargin{ivarargin};
        case 'gainrotationerror'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
        case 'gainvelocityerror'
            ivarargin=ivarargin+1;
            gainVelocityError=varargin{ivarargin};
        case 'flagextendedvector'
            ivarargin=ivarargin+1;
            flagExtendedVector=varargin{ivarargin};
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=disturbance+varargin{ivarargin}(1:12,:);
        case 'wreference'
            ivarargin=ivarargin+1;
            wReference = varargin{ivarargin};
        otherwise
            % Skip data of unused entry
            ivarargin=ivarargin+1;
%             error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[R,w]=rotDyn_stateUnpack(x);
uR=gainRotationError*rot_vee(R,rot_log(R,RReference));
if ~flagVelocityInertiaWeight
    uw= wReference - w;
else
    uw=J*(wReference - w);
end
uw=gainVelocityError*uw;
dwDisturbance=rotDyn_inputUnpack(disturbance);
dw=uR+uw-dwDisturbance;
if flagCancelDynamics
    dw=dw+rotDyn_gyroscopicTerm(w,'inertiaMatrix',J);
end

if flagExtendedVector
    dw=rotDyn_inputPack(dw,'state',R,w);
end

