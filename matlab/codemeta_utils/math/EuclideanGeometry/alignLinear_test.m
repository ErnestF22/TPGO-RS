function alignLinear_test
S=randn(3,10);
M=randn(10,3);
K=randn(3,3);
STarget=K*S;
MTarget=M*K;
[SKEst,STransformed,SCost]=alignLinear(S,STarget);
disp([K SKEst])
disp([STransformed; STarget])
disp(SCost)

[MKEst,MTransformed,MCost]=alignLinear(M,MTarget,'right');
disp([K MKEst])
disp([MTransformed MTarget])
disp(MCost)
