function POCsimilarityTransformLines

x1=real_randGeodFun(randn(2,1),'speed','rand');
x2=real_randGeodFun(randn(2,1));
A=@(t) simTransfEst(x1(t),x2(t));
Ax1=@(t) A(t)*x1(t);
Ax2=@(t) A(t)*x2(t);

figure(1)
subplot(2,1,1)
funPlot(Ax1)
subplot(2,1,2)
funPlot(Ax2)




function A=simTransfEst(x1,x2)
M=[
    x1(1) 0 x1(2) 0 0
    0 x1(1) 0 x1(2) 0
    x2(1) 0 x2(2) 0 -x2(1)
    0 x2(1) 0 x2(2) -x2(2)
];
m=[1;0;0;0];
ak=M\m;

A=reshape(ak(1:4),2,2);