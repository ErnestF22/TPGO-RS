function dw=rotDyn_controlPD(x,varargin)
flagCancelDynamics=false;
flagExtendedVector=false;
RReference=eye(3);
J=eye(3);
gainRotationError=1;
gainVelocityError=1;
flagVelocityInertiaWeight=false;
disturbance=zeros(12,size(x,2));
flagAugmentedSystem = false;

%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
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
        case 'augmentedsystem'
            % controller is following a reference trajectory on
            % TSO(3)xSO(3). RReference should be read from x(13:end)
            flagAugmentedSystem = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

if flagAugmentedSystem
    [R,w,RReference] = rotDyn_stateUnpack(x,'augmentedsystem');
else
    [R,w]=rotDyn_stateUnpack(x);
end

uR=gainRotationError*rot_vee(R,rot_log(R,RReference));
if ~flagVelocityInertiaWeight
    uw=-w;
else
    uw=-J*w;
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
