function R=rotDyn_integrateVelocity(t,w,R0,varargin)
flagRight=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'flagright'
            ivarargin=ivarargin+1;
            flagRight=varargin{ivarargin};
        case 'left'
            flagRight=false;
        case 'right'
            flagRight=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NIt=length(t);
dt=diff(t);
R=zeros(3,3,NIt);
R(:,:,1)=R0;
for it=2:NIt
    RCurrent=R(:,:,it-1);
    if flagRight
        R(:,:,it)=rot_exp(RCurrent,rot_hat(RCurrent,dt(it-1)*w(:,it-1)));
    else
        R(:,:,it)=rot_exp(RCurrent',rot_hat(RCurrent',dt(it-1)*w(:,it-1)))';
    end
end
