function essential_vee_test
Q=essential_randn();
v=essential_randTangentNormVector(Q);
vVec=essential_vee(Q,v);
disp(max(max(abs(v-essential_hat(Q,vVec)))))
Q2=essential_exp(Q,v);
disp(max(max(abs(Q2-[Q(1:3,:)*rot(vVec(1:3)); Q(4:6,:)*rot(vVec(4:6))]))))


