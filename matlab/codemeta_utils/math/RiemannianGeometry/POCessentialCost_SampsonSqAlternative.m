function POCessentialCost_SampsonSqAlternative
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

    funCheckDer(@(t) derDderQAlternative(Q(t),hatv0,x1,x2))

    %Test symmetry
    t=0;
    [~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E(t),x1,x2,'symmetricHess');
    hessOpQ=hessE2hessQ(E(t),gradc,hessOpc);
    
    asym=@(A) (A-A')/2;
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

function [dc,ddc]=derDderQAlternative(Q,hatv0,x1,x2)
    E=essential_toE(Q);
    hatv01=hatv0(1:3,:);
    hatv02=hatv0(4:6,:);
    dE= hatv01'*E+E*hatv02;
    ddE= (hatv01^2)'*E+2*hatv01'*E*hatv02+E*hatv02^2;
    symm=@(A) (A+A')/2;
    
    [~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x1,x2,'symmetricHess');
    dc=trace(gradc'*dE);
    %ddc=trace(dE'*hessOpc(dE))+trace(ddE'*gradc);
    %ddc=trace(dE'*hessOpc(dE))+trace(gradc'*((hatv01^2)'*E+2*hatv01'*E*hatv02+E*hatv02^2));
    %ddc=trace(dE'*hessOpc(dE))+trace(gradc'*(hatv01^2)'*E)+2*trace(gradc'*hatv01'*E*hatv02)+trace(gradc'*E*hatv02^2);
    %ddc=trace(dE'*hessOpc(dE))-trace(hatv01'*symm(E*gradc')*hatv01)+2*trace(gradc'*hatv01'*E*hatv02)-trace(hatv02'*symm(gradc'*E)*hatv02);
    hessEucl=hessOpc(dE);
    %ddc=trace(dE'*hessEucl)-trace(hatv01'*symm(E*gradc')*hatv01)+2*trace(gradc'*hatv01'*E*hatv02)-trace(hatv02'*symm(gradc'*E)*hatv02);
    %ddc=trace((hatv01'*E+E*hatv02)'*hessEucl)-trace(hatv01'*symm(E*gradc')*hatv01)+2*trace(gradc'*hatv01'*E*hatv02)-trace(hatv02'*symm(gradc'*E)*hatv02);
    %ddc=trace(E'*hatv01*hessEucl)+trace(hatv02'*E'*hessEucl)-trace(hatv01'*symm(E*gradc')*hatv01)+2*trace(gradc'*hatv01'*E*hatv02)-trace(hatv02'*symm(gradc'*E)*hatv02);
    %ddc=-trace(hatv01'*hessEucl*E')+trace(hatv02'*E'*hessEucl)-trace(hatv01'*symm(E*gradc')*hatv01)+2*trace(gradc'*hatv01'*E*hatv02)-trace(hatv02'*symm(gradc'*E)*hatv02);
    %ddc=trace(hatv01'*(-hessEucl*E'-symm(E*gradc')*hatv01))+trace(hatv02'*E'*hessEucl)+trace(hatv01'*E*hatv02*gradc')+trace(hatv02'*E'*hatv01*gradc)-trace(hatv02'*symm(gradc'*E)*hatv02);
    %ddc=trace(hatv01'*(-hessEucl*E'-symm(E*gradc')*hatv01+E*hatv02*gradc'))+trace(hatv02'*(E'*hessEucl-symm(gradc'*E)*hatv02+E'*hatv01*gradc));
    ddc=trace(hatv0'*hessOpQBase(hatv0,E,gradc,hessOpc));
end

function hessOpQ=hessE2hessQ(E,gradc,hessOpc)
    hessOpQ=@(hatv0) hessOpQBase(hatv0,E,gradc,hessOpc);
end


function hc=hessOpQBase(hatv0,E,gradc,hessOpc)
    hatv01=hatv0(1:3,:);
    hatv02=hatv0(4:6,:);
    dE= hatv01'*E+E*hatv02;
    hessEucl=hessOpc(dE);
    symm=@(A) (A+A')/2;
    asym=@(A) (A-A')/2;
    hc=[
        -asym(hessEucl*E')-symm(E*gradc')*hatv01+E*hatv02*gradc';
        asym(E'*hessEucl)-symm(gradc'*E)*hatv02+E'*hatv01*gradc
        ];
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


