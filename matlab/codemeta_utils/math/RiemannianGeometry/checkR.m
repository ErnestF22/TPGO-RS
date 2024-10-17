%function b=checkR(R)
%Checks if R is a rotation matrix, b can be true or false
function b=checkR(R)
if(norm(R'*R-eye(3),'fro')>1e-6 || det(R)<0)
    b=false;
else
    b=true;
end

