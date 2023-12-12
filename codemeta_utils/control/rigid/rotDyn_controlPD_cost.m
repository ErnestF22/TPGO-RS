function [phi,gradphi]=rotDyn_controlPD_cost(R,w,varargin)
RReference=eye(3);
J=eye(3);
flagComputeGrad=nargout>1;
gainRotationError=1;
methodGradient='packed';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rreference'
            ivarargin=ivarargin+1;
            RReference=varargin{ivarargin};
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
        case 'gainrotationerror'
            ivarargin=ivarargin+1;
            gainRotationError=varargin{ivarargin};
        case 'methodgradient'
            ivarargin=ivarargin+1;
            methodGradient=lower(varargin{ivarargin});
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NR=size(R,3);
if size(RReference,3)==1 && NR>1
    RReference=repmat(RReference,[1 1 NR]);
end

phi=gainRotationError*rot_dist(R,RReference,'vector').^2/2 ...
    ...% the following is the vectorized version of gainVelocityError*w'*J*w/2
    +squeeze(multiprod(multitransp(permute(w,[1 3 2])),permute(J*w,[1 3 2])))/2;
phi=phi';

if flagComputeGrad
    switch methodGradient
        case 'packed'
            gradphi=rotDyn_statePack(-gainRotationError*rot_log(R,RReference),J*w);
        case 'localcoordinates'
            gradphi=[-gainRotationError*rot_vee(R,rot_log(R,RReference));J*w];
        otherwise
            error('Method for expressing the gradient not recognized.')
    end
end
