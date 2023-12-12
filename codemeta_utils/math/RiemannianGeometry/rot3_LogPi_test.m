function rot3_LogPi_test

plotfuntrials(f,1000)

function c=f()
v=cnormalize(randn(3,1))*pi;
R=rot(v);
c=norm(R-rot(rot3_LogPi(R)),'fro');
