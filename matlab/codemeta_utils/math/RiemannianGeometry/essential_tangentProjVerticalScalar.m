%Gives the scalar value of the projection of a tangent vector on the vertical space
%function p=essential_tangentProjVerticalScalar(Q,v)
function p=essential_tangentProjVerticalScalar(Q,v)
p=0;
for ip=[1 4]
    p=p+squeeze(multiprod(Q(ip+2,:,:),permute(rot_vee(Q(ip:ip+2,:,:),v(ip:ip+2,:,:)),[1 3 2])));
end
