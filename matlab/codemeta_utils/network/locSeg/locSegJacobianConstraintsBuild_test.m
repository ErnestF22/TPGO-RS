function locSegJacobianConstraintsBuild_test
%resetRands()
%constraintType='fixedBearings';
%constraintType='relativeBearings';
%constraintType='distances';
constraintType='relativeRotations';

for dimData=2:3
    [xi,~,~,dxi]=real_randGeodFun(randn(dimData,1));
    [xj,~,~,dxj]=real_randGeodFun(randn(dimData,1));
    [Ri,~,~,~,vi]=rot_randGeodFun(dimData);
    [Rj,~,~,~,vj]=rot_randGeodFun(dimData,'speed',rand);
    Rij=Ri(0)*Rj(0)';

    xixj=@(t) xi(t)-xj(t);
    nxixj=@(t) norm(xi(t)-xj(t));
    tij=@(t) Ri(t)'*xixj(t)/nxixj(t);
    tij0=tij(0);
    dij0=nxixj(0);
    v=[dxi;vi;dxj;vj];

    u=randn(dimData,1);
    switch constraintType
        case 'fixedBearings'
            f=@(t) u'*(xixj(t)/nxixj(t)-tij0);
            Jf=@(t) locSegJacobianConstraintsBuild(constraintType,xi(t),xj(t));
            df=@(t) u'*Jf(t)*v;
        case 'relativeBearings'
            f=@(t) u'*(xixj(t)/nxixj(t)-Ri(t)*tij0);
            Jf=@(t) locSegJacobianConstraintsBuild(constraintType,xi(t),xj(t),Ri(t),tij0);
            df=@(t) u'*Jf(t)*v;
        case 'distances'
            f=@(t) nxixj(t)-dij0;
            Jf=@(t) locSegJacobianConstraintsBuild(constraintType,xi(t),xj(t));
            df=@(t) Jf(t)*v;
        case 'relativeRotations'
            f=@(t) reshape(Ri(t)-Rij*Rj(t),[],1);
            Jf=@(t) locSegJacobianConstraintsBuild(constraintType,xi(t),xj(t),Ri(t),Rj(t),Rij);
            df=@(t) Jf(t)*v;
    end

    figure(dimData-1)
    check_der(f,df)
end