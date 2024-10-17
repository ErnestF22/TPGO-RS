%Compute the projection of a matrix on the tangent space
%function v=essential_tangentProj(Q,v)
function v=essential_tangentProj(Q,v)
v=[rot_tangentProj(Q(1:3,:),v(1:3,:)); rot_tangentProj(Q(4:6,:),v(4:6,:))];
p=essential_tangentProjVerticalScalar(Q,v);
v=v-p/2*[Q(1:3,:)*hat(Q(3,:)');Q(4:6,:)*hat(Q(6,:)')];
