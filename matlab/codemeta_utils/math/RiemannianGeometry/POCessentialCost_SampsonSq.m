function POCessentialCost_SampsonSq
    ez=[0;0;1];
    ezhat=hat(ez);

    [Q,dQ,~,~,v0Vec]=essential_randGeodFun();
    R1=@(t) essential_getR1(Q(t));
    R2=@(t) essential_getR2(Q(t));
    dR1=@(t) essential_getR1(dQ(t));
    dR2=@(t) essential_getR2(dQ(t));
    hatv01=hat(v0Vec(1:3));
    hatv02=hat(v0Vec(4:6));
    hatv0=[hatv01;hatv02];
    ddR1=@(t) dR1(t)*hatv01;
    ddR2=@(t) dR2(t)*hatv02;

    E=@(t) essential_toE(Q(t));
    dE=@(t) hatv01'*E(t)+E(t)*hatv02;
    ddE=@(t) (hatv01^2)'*E(t)+2*hatv01'*E(t)*hatv02+E(t)*hatv02^2;

    %funCheckDer(E,dE)
    %funCheckDer(dE,ddE)

    Nx=3;
    x1=randn(2,Nx);
    x2=randn(2,Nx);
    x=cat(3,x1,x2);

    %funCheckDer(@(t) funDer(E(t),dE(t),x1,x2))
    %funCheckDer(@(t) derDder(E(t),dE(t),ddE(t),x1,x2))

    funCheckDer(@(t) funDerQ(Q(t),hatv0,x))
    funCheckDer(@(t) gradDgradQ(Q(t),hatv0,x))
    funCheckDer(@(t) derDderQ(Q(t),hatv0,x))

    %Test symmetry
    asym=@(A) (A-A')/2;
    [~,~,hessOpQ]=essential_costSampsonSq(Q(0),x,'symmetricHess');
    dQ1=[asym(randn(3)); asym(randn(3))];
    dQ2=[asym(randn(3)); asym(randn(3))];
    disp([trace(dQ1'*hessOpQ(dQ2)) trace(dQ2'*hessOpQ(dQ1))])
end

function [c,dc]=funDerQ(Q,hatv0,x)
    [c,gradcQ]=essential_costSampsonSq(Q,x);
    dc=trace(gradcQ'*hatv0);
end

function [gradcQ,dgradcQ]=gradDgradQ(Q,hatv0,x)
    [~,gradcQ,hessOpcQ]=essential_costSampsonSq(Q,x);
    dgradcQ=hessOpcQ(hatv0);
end

function [dc,ddc]=derDderQ(Q,hatv0,x)
    [~,gradcQ,hessOpcQ]=essential_costSampsonSq(Q,x,'symmetricHess');
    dc=trace(gradcQ'*hatv0);
    ddc=trace(hatv0'*hessOpcQ(hatv0));
end

function [c,dc]=funDer(E,dE,x1,x2)
    [c,gradc]=epipolarCostFromE_SampsonSq(E,x1,x2);
    dc=trace(gradc'*dE);
end

function [dc,ddc]=derDder(E,dE,ddE,x1,x2)
    [~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x1,x2,'symmetricHess');
    dc=trace(gradc'*dE);
    ddc=trace(dE'*hessOpc(dE))+trace(ddE'*gradc);
end

function [gradc,dgradc]=gradDgrad(E,dE,x1,x2)
    [~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x1,x2);
    dgradc=hessOpc(dE);
end


