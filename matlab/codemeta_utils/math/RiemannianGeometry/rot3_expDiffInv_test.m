function rot3_expDiffInv_test

exp=2;
switch exp
    case 1
        plotfuntrials(@ftest,10)
    case 2
        %S=orth(randn(3)); S=S*det(S);
        S=eye(3);
        v=pi*rand*rot_randTangentNormVector(S);
        R0=rot_exp(S,v);
        watR0=rot_randTangentNormVector(R0);
        R=@(t) rot_exp(R0,t*watR0);
        watR=@(t) R(t)*R0'*watR0;

        w=@(t) rot3_expDiffInv(S,R(t),watR(t),'rot');
        
        approx_der(w,0)
end


function r=ftest()
S=orth(randn(3)); S=S*det(S);
v=pi*rand*rot_randTangentNormVector(S);
w=pi*rand*rot_randTangentNormVector(S);

w1=rot3_expDiffInv(S,v,rot3_expDiff(S,v,w));

r=norm(w-w1,'fro');
