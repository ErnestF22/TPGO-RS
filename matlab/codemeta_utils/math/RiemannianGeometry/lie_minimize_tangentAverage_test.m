function lie_minimize_tangentAverage_test
%resetRands();
%testName='costTest';
testName='averagel1';

switch 2
    case 1
        y0=eye(3);
        yi=rot_randn(y0,1,40);
        lf=rot_funs();
    case 2
        y0=essential_eye();
        yi=essential_randn(y0,1,5);
        lf=essential_signed_funs();
end

switch testName
    case 'averagel2'
        fCost=@(y) l2cost(lf,y,yi);
        yInit=yi(:,:,1);
        [y,output]=lie_minimize_tangentAverage(lf,yInit,yi,'displayIt','maxit',25);
        y2=lie_minimizeGradNewton(lf,fCost,yInit,'gradientOnly','displayIt');
        disp([y0 y y2])
        figure(1)
        plot(output.t,output.c)
    case 'averagel1'
        yInit=lie_minimize_tangentAverage_getInitialPoint(lf,yi,'norm',1,'method','twominaverage');
        [y,output]=lie_minimize_tangentAverage(lf,yInit,yi,'displayIt','maxit',25,'norm',1);
        disp([y0 y])
        figure(1)
        plot(output.t,output.c)
    case 'costTest'
        v=lf.randTangentNormVector(y0);
        vVec=lf.vee(y0,v);
        y=@(t) lf.exp(y0,t*v);
        fGeod=@(t) l2costAndDer(lf,y(t),yi,vVec);
        check_der(fGeod,'function',linspace(-1,1,31))
end
        
function [c,gradc]=l2cost(lf,y,yi)
c=0.5*sum(lf.dist(y,yi).^2);
gradc=-sum(lf.vee(y,lf.log(y,yi)),2);

function [c,dc]=l2costAndDer(lf,y,yi,vVec)
[c,gradc]=l2cost(lf,y,yi);
dc=gradc'*vVec;
