%function b=checkg(g)
%Checks if g is a rigid body transformation, b can be true or false
function b=checkg(g)
if(g(4,1)~=0 || g(4,2)~=0 || g(4,3)~=0 || g(4,4)~=1 || ~checkR(g(1:3,1:3)))
    b=false;
else
    b=true;
end
