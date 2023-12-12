function v=quat_rotateVec(q,v)
p=[zeros(1,size(v,2));v];
p=quat_mult(q,quat_mult(p,quat_conj(q)));
v=p(2:4,:);