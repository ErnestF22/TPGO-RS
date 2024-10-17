%generate a set of rotations by chosing a basis in T_e SO(n) and moving of
%pi/2 along each direction in the basis + the identity e is appended to the
%the others
function R=rot_anchors(e)
T=rot_tangentBasis(e);
K=size(T,3);
R=zeros([size(e) K+1]);
for k=1:K
    R(:,:,k)=rot_exp(e,pi*T(:,:,k));
end
R(:,:,K+1)=e;
for k=K+2:K+K+2
    R(:,:,k)=orth(randn(size(e)));
    R(:,:,k)=R(:,:,k)/sign(det(R(:,:,k)));
end
    