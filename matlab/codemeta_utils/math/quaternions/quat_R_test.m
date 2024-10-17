function quat_R_test

plotfuntrials(@fTest)

function e=fTest()
R=rot_randn();
e=norm(R-quat_toR(quat_fromR(R)),'fro');
