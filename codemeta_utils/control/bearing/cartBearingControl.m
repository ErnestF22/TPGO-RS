function u=cartBearingControl(theta,YEval,YGoal,funs,varargin)
methodLinearSpeed='forward-backward';
dmax=10;
knu=0.1;
maxLinearSpeed=Inf;
maxAngularSpeed=Inf;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'dmax'
            ivarargin=ivarargin+1;
            dmax=varargin{ivarargin};
        case 'knu'
            ivarargin=ivarargin+1;
            knu=varargin{ivarargin};
        case 'maxlinearspeed'
            ivarargin=ivarargin+1;
            maxLinearSpeed=varargin{ivarargin};
        case 'maxangularspeed'
            ivarargin=ivarargin+1;
            maxAngularSpeed=varargin{ivarargin};
        case 'methodlinearspeed'
            ivarargin=ivarargin+1;
            methodLinearSpeed=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
dXEval=[cos(theta); sin(theta)];
[bdtheta,g]=bearingCostGeneral_desiredAngle_derBound(YEval,YGoal,funs,dXEval,dmax);

%compute desired linear speed
u=zeros(2,1);
ngSq=g'*g;
switch lower(methodLinearSpeed)
    case 'forward'
        u(1)=sqrt(ngSq);
    case 'forward-backward'
        u(1)=-dXEval'*g;
    otherwise
        error('Cart linear speed method not recognized')
end

%select upper or lower bound according to sign of the angle error
thetad=atan2(-g(2),-g(1));
% display(thetad)
% display(cnormalize(g))


thetae=modAngle(theta-thetad);
if thetae<0
    bdtheta=bdtheta(2);
else
    bdtheta=bdtheta(1);
end

%get control by scaling angular bound according to magnitude of linear speed
u(2)=bdtheta*abs(u(1));
if thetae<0
    u(2)=u(2)+knu;
else
    u(2)=u(2)-knu;
end

% fprintf('derBound %.9f\n',bdtheta);
% fprintf('abs(v) %.9f\n',abs(u(1)));

%rescale computed inputs to satisfy limits
if abs(u(1))>maxLinearSpeed
    u=u/abs(u(1))*maxLinearSpeed;
end
if abs(u(2))>maxAngularSpeed
    u=u/abs(u(2))*maxAngularSpeed;
end
