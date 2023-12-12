function POCessentialDistance
switch 2
    case 1
        %check that tangent to geodesics stays in horizontal space
        Q0=rot_randn([],[],2);
        v0Vec=rot_vee(Q0,rot_randTangentNormVector(Q0));
        e3=[0;0;1];
        u2=[1;1];
        v0Vec=v0Vec-e3*(e3'*v0Vec*u2)*u2'/2;
        v0=rot_hat(Q0,v0Vec);
        [Qt,vt]=rot_geodFun(Q0,v0);
        vtVec=@(t) rot_vee(Qt(t),vt(t));
        m3=@(t) (e3'*vtVec(t)*u2)/2;

        plotfun(m3)
    case 2
        %minimization problem to find log
        e3=[0;0;1];
        Q10=rot_randn([],[],2);
        Q110=Q10(:,:,1);
        Q120=Q10(:,:,2);
        Q20=rot_randn([],[],2);
        Rzt=@(t) rot(t*e3);
        
        Q21t=@(t) Q20(:,:,1)*Rzt(t);
        Q22t=@(t) Q20(:,:,2)*Rzt(t);
        
        ft=@(t) 0.5*(rot_dist(Q110,Q21t(t))^2+rot_dist(Q120,Q22t(t))^2);
        dft=@(t) e3'*(logrot(Q110'*Q21t(t))+logrot(Q120'*Q22t(t)));
        
        check_der(ft,dft,'angle')
    case 3
        %check that tangent to geodesics stays in horizontal space
        Q0=rot_randn([],[],2);
        v0=rot_randTangentNormVector(Q0);

        d=metricInVSpace(Q0,v0);
        
        e3=[0;0;1];
        R01=Q0(:,:,1);
        R02=Q0(:,:,2);
        v0Vec=rot_vee(Q0,v0);
        v0Vec(:,1)=v0Vec(:,1)-d/2*R01'*e3;
        v0Vec(:,2)=v0Vec(:,2)-d/2*R02'*e3;
        v0=rot_hat(Q0,v0Vec);
        
        [Qt,vt]=rot_geodFun(Q0,v0);
        m3t=@(t) metricInVSpace(Qt(t),vt(t));

        m3Logt=@(t) metricInVSpaceLogQ(Q0,Qt(t));
        
        subplot(2,1,1)
        plotfun(m3t)
        subplot(2,1,2)
        plotfun(m3Logt)
    case 4
        %minimization problem to find log
        e3=[0;0;1];
        Q10=rot_randn([],[],2);
        Q110=Q10(:,:,1);
        Q120=Q10(:,:,2);
        Q20=rot_randn([],[],2);
        Rzt=@(t) rot(t*e3);
        
        Q21t=@(t) Rzt(t)*Q20(:,:,1);
        Q22t=@(t) Rzt(t)*Q20(:,:,2);
        
        ft=@(t) 0.5*(rot_dist(Q110,Q21t(t))^2+rot_dist(Q120,Q22t(t))^2);
        dft=@(t) e3'*Q21t(t)*logrot(Q110'*Q21t(t))+e3'*Q22t(t)*logrot(Q120'*Q22t(t));
        dfbt=@(t) e3'*Q110*logrot(Q110'*Q21t(t))+e3'*Q120*logrot(Q120'*Q22t(t));
        ddft=@(t) e3'*Q21t(t)*rot3_expDiffInvMat(eye(3),Q110'*Q21t(t))*Q21t(t)'*e3 ...
            +e3'*Q22t(t)*rot3_expDiffInvMat(eye(3),Q120'*Q22t(t))*Q22t(t)'*e3;
        ddfbt=@(t) e3'*Q110*rot3_expDiffInvMat(eye(3),Q110'*Q21t(t))*Q110'*e3 ...
            +e3'*Q120*rot3_expDiffInvMat(eye(3),Q120'*Q22t(t))*Q120'*e3;
        check_der(ft,dfbt,'angle')
        check_der(dft,ddft,'angle')
        
end

function vVec=Log(Q0,Q1)
vVec=[logrot(Q0(:,:,1)'*Q1(:,:,1)) logrot(Q0(:,:,2)'*Q1(:,:,2))];

function m=metricInVSpaceLogQ(Q0,Q1)
vVec=Log(Q0,Q1);
R1=Q0(:,:,1);
R2=Q0(:,:,2);
m=R1(3,:)*vVec(:,1)+R2(3,:)*vVec(:,2);


function m=metricInVSpace(Q,v)
R1=Q(:,:,1);
R2=Q(:,:,2);
vVec=rot_vee(Q,v);

m=R1(3,:)*vVec(:,1)+R2(3,:)*vVec(:,2);
