function essential_evaluateEpipolarCostSampsonSq_test

    [Q,~,~,~,v0Vec]=essential_randGeodFun();
    hatv01=hat(v0Vec(1:3));
    hatv02=hat(v0Vec(4:6));
    hatv0=[hatv01;hatv02];

    Nx=10;
    x1=randn(2,Nx);
    x2=randn(2,Nx);
    x=cat(3,x1,x2);

    funCheckDer(@(t) funDerQ(Q(t),hatv0,x))
    funCheckDer(@(t) gradDgradQ(Q(t),hatv0,x))
    funCheckDer(@(t) derDderQ(Q(t),hatv0,x))

    %Test symmetry
    [~,~,hessOpQ]=essential_evaluateEpipolarCostSampsonSq(Q(rand),x,'symmetricHess');
    asym=@(A) (A-A')/2;
    dQ1=[asym(randn(3)); asym(randn(3))];
    dQ2=[asym(randn(3)); asym(randn(3))];
    disp([trace(dQ1'*hessOpQ(dQ2)) trace(dQ2'*hessOpQ(dQ1))])
end

function [c,dc]=funDerQ(Q,hatv0,x)
    [c,gradcQ]=essential_evaluateEpipolarCostSampsonSq(Q,x);
    dc=trace(gradcQ'*hatv0);
end

function [gradcQ,dgradcQ]=gradDgradQ(Q,hatv0,x)
    [~,gradcQ,hessOpcQ]=essential_evaluateEpipolarCostSampsonSq(Q,x);
    dgradcQ=hessOpcQ(hatv0);
end

function [dc,ddc]=derDderQ(Q,hatv0,x)
    [~,gradcQ,hessOpcQ]=essential_evaluateEpipolarCostSampsonSq(Q,x,'symmetricHess');
    dc=trace(gradcQ'*hatv0);
    ddc=trace(hatv0'*hessOpcQ(hatv0));
end

function [c,dc]=funDer(E,dE,x1,x2)
    [c,gradc]=essential_evaluateEpipolarCostSampsonSq(E,x1,x2);
    dc=trace(gradc'*dE);
end

function [dc,ddc]=derDder(E,dE,ddE,x1,x2)
    [~,gradc,hessOpc]=essential_evaluateEpipolarCostSampsonSq(E,x1,x2,'symmetricHess');
    dc=trace(gradc'*dE);
    ddc=trace(dE'*hessOpc(dE))+trace(ddE'*gradc);
end

function [gradc,dgradc]=gradDgrad(E,dE,x1,x2)
    [~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x1,x2);
    dgradc=hessOpc(dE);
end


