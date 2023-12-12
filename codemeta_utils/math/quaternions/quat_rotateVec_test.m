function quat_rotateVec_test
v=randn(3,1);
R=rot_randn();

q=quat_fromR(R);

disp([R*v quat_rotateVec(q,v)])
