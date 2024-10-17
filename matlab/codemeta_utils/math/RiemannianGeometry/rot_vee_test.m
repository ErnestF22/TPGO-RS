function rot_vee_test
d=10;
R=rot_randn(eye(d));
dM=d*(d-1)/2;
v=randn(dM,5);
disp(v-rot_vee(R,rot_hat(R,v)))
