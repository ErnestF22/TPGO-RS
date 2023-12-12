function [v,Q2r]=essential_vecLog(Q1,Q2,varargin)
flagSigned=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagsigned'
            ivarargin=ivarargin+1;
            flagSigned=varargin{ivarargin};
        case 'signed'
            flagSigned=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

Q2r=essential_closestRepresentative(Q1,Q2,'flagSigned',flagSigned);
v=[
    logrot(Q1(1:3,:)'*Q2r(1:3,:));
    logrot(Q1(4:6,:)'*Q2r(4:6,:))
  ];

