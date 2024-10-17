function DLog=essential_logDiff(Q1,Q2,varargin)
flagAlreadyClosestRepresentative=false;
components='12';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'closestrepresentative'
            flagAlreadyClosestRepresentative=true;
        case 'components'
            ivarargin=ivarargin+1;
            components=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~flagAlreadyClosestRepresentative
    Q2=essential_closestRepresentative(Q1,Q2);
end

Ra1=essential_getR1(Q1);
Ra2=essential_getR2(Q1);
Rb1=essential_getR1(Q2);
Rb2=essential_getR2(Q2);

Razb1=Ra1'*Rb1;
Razb2=Ra2'*Rb2;

Log1=logrot(Razb1);
Log2=logrot(Razb2);

v1=Rb1(3,:)';
v2=Rb2(3,:)';

D1=rot3_logDiff(eye(3),Log1,'tangentVec');
D2=rot3_logDiff(eye(3),Log2,'tangentVec');

dtOptDen=v1'*D1*v1+v2'*D2*v2;

if any(components=='1')
    A11=Rb1(3,:)*D1/dtOptDen;
    A21=Rb2(3,:)*D2/dtOptDen;
    DLog1=[(-eye(3)+Rb1(3,:)'*A11)*Razb1' Rb1(3,:)'*A21*Razb2';
     Rb2(3,:)'*A11*Razb1' (-eye(3)+Rb2(3,:)'*A21)*Razb2'];
else
    DLog1=[];
end

if any(components=='2')
    A12=Rb1(3,:)*(-D1+hat(Log1))/dtOptDen;
    A22=Rb2(3,:)*(-D2+hat(Log2))/dtOptDen;
    DLog2=[(eye(3)+Rb1(3,:)'*A12) Rb1(3,:)'*A22;
        Rb2(3,:)'*A12 (eye(3)+Rb2(3,:)'*A22)];
else
    DLog2=[];
end

DLog=blkdiag(D1,D2)*[DLog1 DLog2];

% DLog=blkdiag(D1,D2)*...
%     [(-eye(3)+Rb1(3,:)'*A11)*Razb1' Rb1(3,:)'*A21*Razb2'...
%          (eye(3)+Rb1(3,:)'*A12) Rb1(3,:)'*A22;
%      Rb2(3,:)'*A11*Razb1' (-eye(3)+Rb2(3,:)'*A21)*Razb2'...
%         Rb2(3,:)'*A12 (eye(3)+Rb2(3,:)'*A22)];
