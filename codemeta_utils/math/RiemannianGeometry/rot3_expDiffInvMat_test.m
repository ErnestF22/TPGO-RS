function rot3_expDiffInvMat_test
%resetRands();

expNum=3;
switch expNum
    case 1
        %compare various methods to compute the differential (3 ways for
        %the matrix representation and one for the mapping representation)
        R1=rot_randn();
        R2=rot_randn();
        S=rot_randn();

        v2=rand*rot_randTangentNormVector(R2);
        D1=rot3_expDiffInvMat(R1,R2,'rot','method','inverseClosedForm');
        D2=rot3_expDiffInvMat(R1,R2,'rot','method','inverseSeries');
        D3=rot3_expDiffInvMat(R1,R2,'rot','method','closedForm');

        disp('From inverse/inverse series/closed form')
        disp([D1 D2 D3])

        D=D3;
        disp('Map/Matrix')
        disp([rot3_expDiffInv(R1,R2,v2,'rot') rot_hat(R1,D*rot_vee(R2,v2))])
        disp('Map/Matrix with conjugate through S')
        disp([rot3_expDiffInv(R1,R2,R2*S*R2'*v2*S','rot') rot_hat(R1,D*S*rot_vee(R2,v2))])
    case 2
        %check the definition of differential of log
        R0=rot_randn(eye(3));
        [Rt,vt]=rot_randGeodFun(R0);
        
        LogRt=@(t) logrot(Rt(t));
        thetat=@(t) norm(LogRt(t));
        vVect=@(t) rot_vee(Rt(t),vt(t));
        DLogRt=@(t) rot3_expDiffInvMat(eye(3),Rt(t),'rot');
        dLogRt=@(t) DLogRt(t)*vVect(t);
        
        eigDLogRt=@(t) eig(DLogRt(t));
        eigDLogRtTheory=@(t) [1;thetat(t)*0.5*(1j+cot(thetat(t)/2));thetat(t)*0.5*(-1j+cot(thetat(t)/2))];
        
        figure(1)
        check_der(LogRt,dLogRt,'angle')

        %check eigenvalues
        figure(2)
        subplot(2,1,1)
        plotfun(@(t) real(eigDLogRt(t)),'angle','r')
        hold on
        plotfun(@(t) real(eigDLogRtTheory(t)),'angle','gx')
        hold off
        title('Real part')
        subplot(2,1,2)
        plotfun(@(t) imag(eigDLogRt(t)),'angle','r')
        hold on
        plotfun(@(t) imag(eigDLogRtTheory(t)),'angle','gx')
        hold off
        title('Imaginary part')
        legend('Computed','Expected from theory')
    case 3
        %check property for multiplication by R'
        R=rot_randn();
        DLogR=rot3_expDiffInvMat(eye(3),R,'rot');
        disp('[DLogR*R'' R''*DLogR DLogR'']')
        disp([DLogR*R' R'*DLogR DLogR'])
        disp(['Max abs diff=' num2str(max(max(abs(DLogR*R'-R'*DLogR))))])
        disp(['Max abs diff=' num2str(max(max(abs(DLogR*R'-DLogR'))))])
end
        