function POCCVXLowRank
w=0.1*randn(3,1);
v=randn(3,1);
n=cnormalize(randn(3,1));
MTruth=v*n';
MSymTruth=multisym(MTruth);
H=hat3(w)+MTruth;

w2=0.1*randn(3,1);
v2=randn(3,1);
MTruth2=v2*n';
H2=hat3(w2)+MTruth2;

numExperiment=5;
switch numExperiment
    case 1
        cvx_begin
            variable M(3,3)
            variable W(3,3)
            minimize norm_nuc(M)+0.01*trace(W(:)'*W(:))
            subject to
                M+W==H
                W==-W'
        cvx_end
        disp([M MTruth])
    case 2
        cvx_begin
            variable MSym(3,3)
            minimize norm_nuc(MSym)
            subject to
                (H+H')/2==MSym;
        cvx_end
        disp([MSymTruth MSym])
    case 3
        cvx_begin
            variable MBigSym(6,6) symmetric semidefinite
            minimize trace(MBigSym)
            subject to
                H+H'==MBigSym(1:3,4:6)+MBigSym(4:6,1:3);
                trace(MBigSym(4:6,4:6))==1;
        cvx_end
        MBigSym=full(MBigSym);
        disp([H+H' MBigSym(1:3,4:6)+MBigSym(4:6,1:3)])
        svd(MBigSym)
    case 4
        %order in MBigSym: [v_(1:3);v2_(4:6);n_(7:9)]
        cvx_begin
            variable MBigSym(9,9) symmetric semidefinite
            minimize trace(MBigSym)
            subject to
                H+H'==MBigSym(1:3,7:9)+MBigSym(7:9,1:3);
                H2+H2'==MBigSym(4:6,7:9)+MBigSym(4:6,7:9);
                trace(MBigSym(7:9,7:9))==1;
        cvx_end
        MBigSym=full(MBigSym);
        %disp([H+H' MBigSym(1:3,7:6)+MBigSym(4:6,1:3)])
        disp(H+H')
        disp(H2+H2')
        disp(MBigSym+MBigSym')
        svd(MBigSym)
    case 5
        cvx_begin
            variable MOuter(3,3)
            minimize norm_nuc(MOuter)
            subject to
                H+H'==MOuter(1:3,1:3)+MOuter(1:3,1:3)';
        cvx_end
        disp(H+H')
        disp(MOuter(1:3,1:3)+MOuter(1:3,1:3)')
    case 6
        cvx_begin
            variable MOuter(6,3)
            minimize norm_nuc(MOuter)
            subject to
                H+H'==MOuter(1:3,1:3)+MOuter(1:3,1:3)';
                H2+H2'==MOuter(4:6,1:3)+MOuter(4:6,1:3)';
        cvx_end
        disp(H+H')
        disp(H2+H2')
        disp(MOuter(1:3,1:3)+MOuter(1:3,1:3)')
        disp(MOuter(4:6,1:3)+MOuter(4:6,1:3)')
        [U,S,V]=svd(MOuter);
        disp([[v;v2] S(1,1)*U(:,1)])
        disp([n V(:,1)])
end
keyboard
