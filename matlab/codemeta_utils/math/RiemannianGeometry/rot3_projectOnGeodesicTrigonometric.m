%Project a rotation on a closed geodesic (trigonometric method)
%function [RProj,tProj]=rot3_projectOnGeodesicTrigonometric(R,Rz0,vRz0)
%Inputs
%   R   point to project
%   R0  starting point of the geodesic
%   vR0 tangent vector of the geodesic at R0
%With respect to purely algebraic version, it returns also the parameter
%tProj for the geodesic.
function [RProj,tProj]=rot3_projectOnGeodesicTrigonometric(R,Rz0,vRz0)
if ~exist('Rz0') || isempty(Rz0)
    Rz0=eye(3);
end
if ~exist('vRz0') || isempty(vRz0)
    vRz0=rot_hat(Rz0,[0;0;1]);
end

uz=rot_vee(Rz0,vRz0);
[RCanonical,Re]=rot3_projectOnGeodesic_canonicalTo(R,Rz0,uz);
c1=RCanonical(1,1)+RCanonical(2,2);
c2=RCanonical(1,2)-RCanonical(2,1);

phi=atan2(c1,c2);
tProj=modAngle(phi-pi/2);
RProj=rot_exp(Rz0,tProj*vRz0);
