function vVec=essential_vee(Q,v)
Nv=size(v,3);
vVec=zeros(6,Nv);
for iv=1:Nv
    vVec(1:3,iv)=vee(Q(1:3,:)'*v(1:3,:,iv));
    vVec(4:6,iv)=vee(Q(4:6,:)'*v(4:6,:,iv));
end