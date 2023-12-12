function [Qt,vt]=essential_geodFun(Q0,v0,varargin)
s=1;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'speed'
            ivarargin=ivarargin+1;
            s=varargin{ivarargin};
        case 'randspeed'
            s=rand;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

v0Eye=s*[Q0(1:3,:)'*v0(1:3,:);Q0(4:6,:)'*v0(4:6,:)];
Qt=@(t) essential_exp(Q0,t*v0);
vt=@(t) [essential_getR1(Qt(t))*v0Eye(1:3,:); essential_getR2(Qt(t))*v0Eye(4:6,:)];
