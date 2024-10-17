function POCEpipolarResidualsSampsonSqDerAndHess
x1=rand(2,1);
x2=rand(2,1);

x1=homogeneous(x1,3);
x2=homogeneous(x2,3);


[E,dE,~,~,ddE]=real_randGeodFun(eye(3),'speed','quadratic');

%funCheckDer(@(t) funDer(E(t),dE(t),x1,x2),'function')
%funCheckDer(@(t) derDder(E(t),dE(t),ddE(t),x1,x2),'function')

[E,dE]=real_randGeodFun(eye(3));
funCheckDer(@(t) gradDgrad(E(t),dE(t),x1,x2))


function [grade,dgrade]=gradDgrad(E,dE,x1,x2)
%Residuals (arguments of the squares)
numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

%Numerator and denominator of the Sampson cost
numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;

gradNumeratorRes=x1*x2';

z=zeros(3,1);
graddenominator11Res=[x2 z z]';
graddenominator12Res=[z x2 z]';
graddenominator21Res=[x1 z z];
graddenominator22Res=[z x1 z];

%Gradient of the numerator and the denominator
gradNumerator=2*numeratorRes*gradNumeratorRes;
%gradDenominator=2*([x2*denominator1Res',z]'+[x1*denominator2Res',z]);
gradDenominator=2*(denominator1Res(1)*graddenominator11Res...
    +denominator1Res(2)*graddenominator12Res...
    +denominator2Res(1)*graddenominator21Res...
    +denominator2Res(2)*graddenominator22Res...
    );

%Gradient of the Sampson cost
gradeNumerator=denominator*gradNumerator-numerator*gradDenominator;
gradeDenominator=(denominator^2);
grade=gradeNumerator/gradeDenominator;

%-------------------------------------------------------------------------%

%dnumerator=trace(gradNumerator'*dE);
%ddenominator=trace(gradDenominator'*dE);

%dnumeratorRes=trace(gradNumeratorRes'*dE);

dgradNumerator=2*trace(gradNumeratorRes'*dE)*gradNumeratorRes;
%dGradNumeratorRes is zero, so this term cancels out: +numeratorRes*dgradNumeratorRes

dgradDenominator=2*(trace(dE'*graddenominator11Res)*graddenominator11Res...
    +trace(dE'*graddenominator12Res)*graddenominator12Res...
    +trace(dE'*graddenominator21Res)*graddenominator21Res...
    +trace(dE'*graddenominator22Res)*graddenominator22Res...
    );

dgradeNumerator=trace(gradDenominator'*dE)*gradNumerator+denominator*dgradNumerator...
    -trace(gradNumerator'*dE)*gradDenominator-numerator*dgradDenominator;
dgradeDenominator=2*denominator*trace(gradDenominator'*dE);

dgrade=(dgradeNumerator*gradeDenominator-gradeNumerator*dgradeDenominator)/(gradeDenominator^2);

function [grade,dgrade]=gradDgradRef(E,dE,x1,x2)
%Residuals (arguments of the squares)
numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

%Numerator and denominator of the Sampson cost
numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;

gradNumeratorRes=x1*x2';

z=zeros(3,1);
graddenominator11Res=[x2 z z]';
graddenominator12Res=[z x2 z]';
graddenominator21Res=[x1 z z];
graddenominator22Res=[z x1 z];

%Gradient of the numerator and the denominator
gradNumerator=2*numeratorRes*gradNumeratorRes;
%gradDenominator=2*([x2*denominator1Res',z]'+[x1*denominator2Res',z]);
gradDenominator=2*(denominator1Res(1)*graddenominator11Res...
    +denominator1Res(2)*graddenominator12Res...
    +denominator2Res(1)*graddenominator21Res...
    +denominator2Res(2)*graddenominator22Res...
    );

%Gradient of the Sampson cost
gradeNumerator=denominator*gradNumerator-numerator*gradDenominator;
gradeDenominator=(denominator^2);
grade=gradeNumerator/gradeDenominator;

%-------------------------------------------------------------------------%

dnumerator=trace(gradNumerator'*dE);
ddenominator=trace(gradDenominator'*dE);

dnumeratorRes=trace(gradNumeratorRes'*dE);

dgradNumerator=2*dnumeratorRes*gradNumeratorRes;
%dGradNumeratorRes is zero, so this term cancels out: +numeratorRes*dgradNumeratorRes

dgradDenominator=2*(trace(dE'*graddenominator11Res)*graddenominator11Res...
    +trace(dE'*graddenominator12Res)*graddenominator12Res...
    +trace(dE'*graddenominator21Res)*graddenominator21Res...
    +trace(dE'*graddenominator22Res)*graddenominator22Res...
    );

dgradeNumerator=ddenominator*gradNumerator+denominator*dgradNumerator...
    -dnumerator*gradDenominator-numerator*dgradDenominator;
dgradeDenominator=2*denominator*ddenominator;

dgrade=(dgradeNumerator*gradeDenominator-gradeNumerator*dgradeDenominator)/(gradeDenominator^2);


function [de,dde]=derDder(E,dE,ddE,x1,x2)
%Residuals (arguments of the squares)
numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

%Numerator and denominator of the Sampson cost
numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;

gradNumeratorRes=x1*x2';

z=zeros(3,1);
graddenominator11Res=[x2 z z]';
graddenominator12Res=[z x2 z]';
graddenominator21Res=[x1 z z];
graddenominator22Res=[z x1 z];

%Gradient of the numerator and the denominator
gradNumerator=2*numeratorRes*gradNumeratorRes;
gradDenominator=2*([x2*denominator1Res',z]'+[x1*denominator2Res',z]);

%Gradient of the Sampson cost
gradeNumerator=denominator*gradNumerator-numerator*gradDenominator;
grade=gradeNumerator/(denominator^2);

de=trace(grade'*dE);

%-------------------------------------------------------------------------%

%Numerator and denominator of the derivative of the Sampson cost
deNumerator=trace(gradeNumerator'*dE);
deDenominator=(denominator.^2);

%Second derivative of the numerator and denominator of the Sampson cost
hessOpNumerator=2*trace(gradNumeratorRes'*dE)*gradNumeratorRes;

hessOpDenominator=2*(trace(dE'*graddenominator11Res)*graddenominator11Res...
    +trace(dE'*graddenominator12Res)*graddenominator12Res...
    +trace(dE'*graddenominator21Res)*graddenominator21Res...
    +trace(dE'*graddenominator22Res)*graddenominator22Res...
    );

%Derivative of the numerator and denominator of the derivative of the
%Sampson cost
hessOpdeNumerator=denominator*hessOpNumerator-numerator*hessOpDenominator;

graddeDenominator=2*denominator*gradDenominator;

hessOpde=(deDenominator*hessOpdeNumerator-deNumerator*graddeDenominator)/(deDenominator^2);
dde=trace(hessOpde'*dE)+trace(grade'*ddE);

function [de,dde]=derDderLong(E,dE,ddE,x1,x2)
%Residuals (arguments of the squares)
numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

%Numerator and denominator of the Sampson cost
numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;

gradNumeratorRes=x1*x2';
dnumeratorRes=trace(gradNumeratorRes'*dE);
%gradDenominator1Res=
z=zeros(3,1);
graddenominator11Res=[x2 z z]';
graddenominator12Res=[z x2 z]';
graddenominator21Res=[x1 z z];
graddenominator22Res=[z x1 z];
ddenominator11Res=trace(dE'*graddenominator11Res);
ddenominator12Res=trace(dE'*graddenominator12Res);
ddenominator21Res=trace(dE'*graddenominator21Res);
ddenominator22Res=trace(dE'*graddenominator22Res);

dnumerator=2*numeratorRes.*dnumeratorRes;
ddenominator1=2*trace([x2*denominator1Res',zeros(3,1)]*dE);
ddenominator2=2*trace([x1*denominator2Res',zeros(3,1)]'*dE);
ddenominator=ddenominator1+ddenominator2;

%Gradient of the numerator and the denominator
gradNumerator=2*numeratorRes*gradNumeratorRes;
%dnumerator=trace(gradNumerator'*dE);
gradDenominator=2*([x2*denominator1Res',zeros(3,1)]'+[x1*denominator2Res',zeros(3,1)]);
%ddenominator=trace(gradDenominator'*dE);

%Gradient of the Sampson cost
gradeNumerator=denominator*gradNumerator-numerator*gradDenominator;
grade=gradeNumerator/(denominator^2);

de=trace(grade'*dE);

%-------------------------------------------------------------------------%

%Numerator and denominator of the derivative of the Sampson cost
deNumerator=dnumerator.*denominator-numerator*ddenominator;
deDenominator=(denominator.^2);

%-------------------------------------------------------------------------%

%Second derivative of the residuals
%these all depend only on ddE
ddnumeratorRes=x1'*ddE*x2;
dddenominator1Res=ddE(1:2,:)*x2;
dddenominator2Res=ddE(:,1:2)'*x1;

%Second derivative of the numerator and denominator of the Sampson cost
%ddnumerator=2*trace(gradNumeratorRes'*dE)^2+2*numeratorRes*ddnumeratorRes;
hessOpNumerator=2*trace(gradNumeratorRes'*dE)*gradNumeratorRes;
ddnumerator=trace(hessOpNumerator'*dE)+trace(gradNumerator'*ddE);

%dddenominator1=2*trace(dE'*graddenominator11Res)^2+2*trace(dE'*graddenominator12Res)^2;
%dddenominator2=2*trace(dE'*graddenominator21Res)^2+2*trace(dE'*graddenominator22Res)^2;
hessOpDenominator=2*(trace(dE'*graddenominator11Res)*graddenominator11Res...
    +trace(dE'*graddenominator12Res)*graddenominator12Res...
    +trace(dE'*graddenominator21Res)*graddenominator21Res...
    +trace(dE'*graddenominator22Res)*graddenominator22Res...
    );
dddenominator=trace(hessOpDenominator'*dE)+trace(gradDenominator'*ddE);
% dddenominator1=2*(trace(dE'*graddenominator11Res)^2+trace(dE'*graddenominator12Res)^2+denominator1Res'*dddenominator1Res);
% dddenominator2=2*(trace(dE'*graddenominator21Res)^2+trace(dE'*graddenominator22Res)^2+denominator2Res'*dddenominator2Res);
% dddenominator=dddenominator1+dddenominator2;

%Derivative of the numerator and denominator of the derivative of the
%Sampson cost
hessOpdeNumerator=denominator*hessOpNumerator-numerator*hessOpDenominator;
ddeNumerator=trace(hessOpdeNumerator'*dE)+trace(gradeNumerator'*ddE);
%ddeNumerator=ddnumerator*denominator-numerator*dddenominator; 

graddeDenominator=2*denominator*gradDenominator;
ddeDenominator=trace(graddeDenominator'*dE);
%ddeDenominator=2*denominator*ddenominator;

hessOpde=(deDenominator*hessOpdeNumerator-deNumerator*graddeDenominator)/(deDenominator^2);
dde=trace(hessOpde'*dE)+trace(grade'*ddE);
%dde=(ddeNumerator*deDenominator-deNumerator*ddeDenominator)/(deDenominator^2);

function [de,dde]=derDderRef(E,dE,ddE,x1,x2)
%Residuals (arguments of the squares)
numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

%Numerator and denominator of the Sampson cost
numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;

%Derivative of the residuals
dnumeratorRes=sum(x1.*(dE*x2));
ddenominator1Res=dE(1:2,:)*x2;
ddenominator2Res=dE(:,1:2)'*x1;

%Derivatives of the numerator and denominator of the Sampson cost
dnumerator=2*numeratorRes.*dnumeratorRes;
ddenominator1=2*denominator1Res.*ddenominator1Res;
ddenominator2=2*denominator2Res.*ddenominator2Res;
ddenominator=sum(ddenominator1+ddenominator2);

%Numerator and denominator of the derivative of the Sampson cost
deNumerator=dnumerator.*denominator-numerator*ddenominator;
deDenominator=(denominator.^2);
de=deNumerator/deDenominator;

%Second derivative of the residuals
ddnumeratorRes=x1'*ddE*x2;
dddenominator1Res=ddE(1:2,:)*x2;
dddenominator2Res=ddE(:,1:2)'*x1;

%Second derivative of the numerator and denominator of the Sampson cost
ddnumerator=2*(dnumeratorRes^2+numeratorRes*ddnumeratorRes);
dddenominator1=2*(ddenominator1Res.^2+denominator1Res.*dddenominator1Res);
dddenominator2=2*(ddenominator2Res.^2+denominator2Res.*dddenominator2Res);
dddenominator=sum(dddenominator1+dddenominator2);

%Derivative of the numerator and denominator of the derivative of the
%Sampson cost
ddeNumerator=ddnumerator*denominator-numerator*dddenominator; 
%This cancels out: +dnumerator*ddenominator-dnumerator*ddenominator
ddeDenominator=2*denominator*ddenominator;
dde=(ddeNumerator*deDenominator-deNumerator*ddeDenominator)/(deDenominator^2);



function [e,de]=funDer(E,dE,x1,x2)

numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;
e=numerator./denominator;

%dnumeratorRes=trace(x2*x1'*dE);
%ddenominator1Res=dE(1:2,:)*x2;
%ddenominator2Res=dE(:,1:2)'*x1;

%dnumerator=2*numeratorRes.*dnumeratorRes;
%ddenominator1=2*trace([x2*denominator1Res',zeros(3,1)]*dE);
%ddenominator2=2*trace([x1*denominator2Res',zeros(3,1)]'*dE);
%ddenominator=ddenominator1+ddenominator2;

gradNumerator=2*numeratorRes*x1*x2';
%dnumerator=trace(gradNumerator'*dE);
gradDenominator=2*([x2*denominator1Res',zeros(3,1)]'+[x1*denominator2Res',zeros(3,1)]);
%ddenominator=trace(gradDenominator'*dE);
grade=(denominator*gradNumerator-numerator*gradDenominator)/(denominator^2);


%de=(dnumerator.*denominator-numerator*ddenominator)./(denominator.^2);
de=trace(grade'*dE);

function [e,de]=funDerRef(E,dE,x1,x2)

numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;
e=numerator./denominator;

dnumeratorRes=sum(x1.*(dE*x2));
ddenominator1Res=dE(1:2,:)*x2;
ddenominator2Res=dE(:,1:2)'*x1;

dnumerator=2*numeratorRes.*dnumeratorRes;
ddenominator1=2*denominator1Res.*ddenominator1Res;
ddenominator2=2*denominator2Res.*ddenominator2Res;
ddenominator=sum(ddenominator1+ddenominator2);

de=(dnumerator.*denominator-numerator*ddenominator)./(denominator.^2);


