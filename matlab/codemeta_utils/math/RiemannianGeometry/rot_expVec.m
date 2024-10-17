%function vVec=rot_expVec(R1,v1)
%If the function is called with a single argument, R1 is assumed to be
%equal to the identity
function R2=rot_expVec(R1,vVec)
flagProvidedRotation=nargin>1;

if ~flagProvidedRotation
    vVec=R1;
    R1=[];
end
if isempty(R1)
    R1=eye(3);
end

R2=rot_exp(R1,rot_hat(R1,vVec));

