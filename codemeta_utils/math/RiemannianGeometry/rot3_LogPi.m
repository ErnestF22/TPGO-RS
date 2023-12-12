%Computes the Log of the rotation for the particular case where d(I,R)=pi
function v=rot3_LogPi(R)
thresholdZero=1e-7;

%get the rank-1 component (u*u')
uut=(R+eye(3))/2;
%get magnitudes
vAbs=pi*sqrt(max(diag(uut),0));
q(1)=uut(1,2);
q(2)=uut(1,3);
q(3)=uut(2,3);

%vZeros=num2str(vAbs'<thresholdZero,'%d');
vZeros=char((vAbs'<thresholdZero)+'0');

switch vZeros
    case '011'
        %disp('Rotation about x')
        v=[pi;0;0];
    case '101'
        %disp('Rotation about y')
        v=[0;pi;0];
    case '110'
        %disp('Rotation about z')
        v=[0;0;pi];
    case '100'
        %disp('u(1) is zero')
        if q(3)>0
            v=vAbs;
        else
            v=diag([0 1 -1])*vAbs;
        end
    case '010'
        %disp('u(2) is zero')
        if q(2)>0
            v=vAbs;
        else
            v=diag([1 0 -1])*vAbs;
        end
    case '001'
        %disp('u(3) is zero')
        if q(1)>0
            v=vAbs;
        else
            v=diag([1 -1 0])*vAbs;
        end
    case '000'
        s(1)=1;
        s(2)=sign(q(1))*s(1);
        s(3)=sign(q(2))*s(1);
        v=diag(s)*vAbs;
    otherwise
        error('%s: this case is impossible. Check thresholdZero.',vZeros)
end

