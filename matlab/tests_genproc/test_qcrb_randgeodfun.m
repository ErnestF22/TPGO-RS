function test_qcrb_randgeodfun

qc0=make_rand_stiefel_3d_array(4,4,1);
if det(qc0) < 0
    qc0(:,1) = - qc0(:,1);
end
rb0=make_rand_stiefel_3d_array(2,2,1);
if det(rb0) < 0
    rb0(:,1) = - rb0(:,1);
end

disp("check_is_rotation(qc0)")
disp(check_is_rotation(qc0))
disp("check_is_rotation(rb0)")
disp(check_is_rotation(rb0))


M.qc = qc0;
M.rb = rb0;

[st,dst,s0,ds0,vVec,ddst,dvVec] = qcrb_randGeodFun(M);

disp('st(0).qc')
disp(st(0).qc)
disp('st(0).rb')
disp(st(0).rb)

disp('dst(0).qc')
disp(dst(0).qc)
disp('dst(0).rb')
disp(dst(0).rb)

disp('ddst(0).qc')
disp(ddst(0).qc)
disp('ddst(0).rb')
disp(ddst(0).rb)

end %file function