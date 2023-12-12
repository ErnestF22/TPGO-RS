function rigidDyn_statePackRT_test
NPoses=5;
R=rot_randn(eye(3),[],NPoses);
w=randn(3,NPoses);
T=randn(3,NPoses);
v=randn(3,NPoses);

[R1,w1,T1,v1]=rigidDyn_stateUnpackRT(rigidDyn_statePackRT(R,w,T,v));

disp(norm(R(:)-R1(:)))
disp(norm(w(:)-w1(:)))
disp(norm(T(:)-T1(:)))
disp(norm(v(:)-v1(:)))

