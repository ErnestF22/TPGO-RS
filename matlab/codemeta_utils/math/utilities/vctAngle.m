%function a=vctAngle(x1,x2)
%Returns the angle between x1 and x2 computed using the QR decomposition
function a=vctAngle(x1,x2)

if size(x1,1)==1
    a=pi*(x1==-x2);
else
    [Q,R]=qr([x1 x2]);
    a=atan2(abs(R(2,2)),sign(R(1,1))*R(1,2));
end
