%Builds the essential matrix representation E from its QREM representation
%function E=essential_getE(Q)
function E=essential_getE(Q)
e3hat=[0 -1 0; 1 0 0; 0 0 0];
N=size(Q,3);
E=zeros(3,3,N);
for iN=1:N
    E(:,:,iN)=Q(1:3,:,iN)'*e3hat*Q(4:6,:,iN);
end
