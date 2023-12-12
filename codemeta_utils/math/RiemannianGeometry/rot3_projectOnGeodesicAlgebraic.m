%Project a rotation on a closed geodesic (purely algebraic method)
%function RProj=rot3_projectOnGeodesicAlgebraic(R,Rz0,vRz0)
%Inputs
%   R   point(s) to project
%   R0  starting point of the geodesic
%   vR0 tangent vector of the geodesic at R0
%
function RProj=rot3_projectOnGeodesicAlgebraic(R,Rz0,vRz0)
if ~exist('Rz0','var') || isempty(Rz0)
    Rz0=eye(3);
end
if ~exist('vRz0','var') || isempty(vRz0)
    vRz0=rot_hat(Rz0,[0;0;1]);
end

NR=size(R,3);
if NR>1
    RProj=zeros(size(R));
    for iR=1:NR
        RProj(:,:,iR)=rot3_projectOnGeodesicAlgebraic(R(:,:,iR),Rz0,vRz0);
    end
else
    uz=rot_vee(Rz0,vRz0);
    %[RCanonical,Re]=rot3_projectOnGeodesic_canonicalTo(R,Rz0,uz);
    %RProjCanonical=householderRotation3Min(RCanonical(:,3),3)'*RCanonical;
    %RProj=rot3_projectOnGeodesic_canonicalFrom(RProjCanonical,Rz0,Re);
    H=householderRotation(uz,R'*Rz0*uz);
    Ruzpi=2*(uz*uz')-eye(3);
    RProj=R*H*Ruzpi;
end
