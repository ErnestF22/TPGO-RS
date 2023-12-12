function POCbearingCostRotationVisibility
%resetRands(1)
y0=[1;0];
NLandmarks=5;
offset=5*[1;1];
xLandmarks=randn(2,NLandmarks)+offset*ones(1,NLandmarks);
x0=[0;0];
%rMax=12.5;

funs=bearingCostFunctions('cosine');

S=[0 -1; 1 0];
R0=eye(2);
dRVec=rand;

R=@(t) rot_exp(R0,t*dRVec*S);
[T,~,~,dT]=real_randGeodFun(x0,'speed',rand);
t=linspace(-pi,pi);
ny=@(t) bearingComputeRanges(T(t),xLandmarks);
rMax=max(max(evalfun(ny,t)));

figure(1)
subplot(3,1,1)
check_der(@(t) costAndDer(R(t),dRVec,T(t),dT,xLandmarks,y0,funs),'function',t)

dPhiBound=@(t) derBound(R(t),dRVec,T(t),dT,xLandmarks,y0,funs,rMax);
dPhiBoundR=@(t) gradRBound(R(t),T(t),xLandmarks,y0,funs,rMax)*dRVec;
dPhiControl=@(t) derAlongControl(R(t),T(t),xLandmarks,y0,funs,rMax);
hold on
plotfun(dPhiBound,t,'c')
plotfun(dPhiBoundR,t,'m')
hold off
subplot(3,1,2)
plotfun(ny,t);
hold on
plotfun(@(t) rMax,t,'r:');
hold off
subplot(3,1,3)
plotfun(dPhiControl,t)

function [phi,dPhi]=costAndDer(R,dRVec,x,dx,xLandmarks,y0,funs)
NLandmarks=size(xLandmarks,2);
[y,ny]=bearingCompute(x,xLandmarks);
S=[0 -1; 1 0];
c=bearingComputeCosine(R*y,y0);
phi=ny*funs.f(c)';
gradRPhiVec=0;
for iLandmark=1:NLandmarks
    gradRPhiVec=gradRPhiVec+ny(iLandmark)*funs.df(c(iLandmark))*y0'*R*S*y(:,iLandmark);
end
gradxPhi=bearingCostGeneral_gradient(y,R'*y0,funs);

dPhi=gradRPhiVec*dRVec+gradxPhi'*dx;

function d=derAlongControl(R,x,xLandmarks,y0,funs,rMax)
[dRVec,dx]=control(R,x,xLandmarks,y0,funs,rMax);
[~,d]=costAndDer(R,dRVec,x,dx,xLandmarks,y0,funs);

function [dRVec,dx]=control(R,x,xLandmarks,y0,funs,rMax)
dMax=10;
gradRPhiBound=gradRBound(R,x,xLandmarks,y0,funs,rMax);
flagGradRSignAvailable=false;
if all(gradRPhiBound>=0)
    %gradRPhi is positive
    b=gradRPhiBound(1);
    flagGradRSignAvailable=true;
end
if all(gradRPhiBound<=0)
    %gradRPhi is negative
    b=gradRPhiBound(2);
    flagGradRSignAvailable=true;
end

if flagGradRSignAvailable
    dRVec=-dMax*b/b^2;
    dx=0;
else
    dRVec=0;
    y=bearingCompute(x,xLandmarks);
    g=bearingCostGeneral_gradient(y,R'*y0,funs);
    dx=-dMax/(g'*g)*g;
end

function dPhiBound=derBound(R,dRVec,x,dx,xLandmarks,y0,funs,rMax)
flagComputeExactDer=false;
NLandmarks=size(xLandmarks,2);
y=bearingCompute(x,xLandmarks);

gradxPhi=bearingCostGeneral_gradient(y,R'*y0,funs);
ax=gradxPhi'*dx;
gradRPhiBound=gradRBound(R,x,xLandmarks,y0,funs,rMax)*dRVec;
dPhiBound=gradRPhiBound+ax;

if flagComputeExactDer
    ny=bearingComputeRanges(x,xLandmarks);
    S=[0 -1; 1 0];
    c=bearingComputeCosine(R*y,y0);
    ab=zeros(1,NLandmarks);
    for iLandmark=1:NLandmarks
        ab(iLandmark)=ny(iLandmark)*funs.df(c(iLandmark))*y0'*R*S*y(:,iLandmark);
    end
    ab=ab*dRVec;
    dPhi=sum(ab)+ax;
end

%Compute upper and lower bounds on the gradient w.r.t. R
%Outputs
%   gradRPhiBound(1)   upper bound
%   gradRPhiBound(2)   lower bound
function gradRPhiBound=gradRBound(R,x,xLandmarks,y0,funs,rMax)
NLandmarks=size(xLandmarks,2);
[y,ny]=bearingCompute(x,xLandmarks);
S=[0 -1; 1 0];
c=bearingComputeCosine(R*y,y0);
a=zeros(1,NLandmarks);
for iLandmark=1:NLandmarks
    a(iLandmark)=rMax*funs.df(c(iLandmark))*y0'*R*S*y(:,iLandmark);
end
%a=a*dRVec;
gradRPhiBound(1)=zeroIfEmpty(sum(a(a>0)));
gradRPhiBound(2)=zeroIfEmpty(sum(a(a<0)));

%Return zero if a is empty, and a itself otherwise
function a=zeroIfEmpty(a)
if isempty(a)
    a=0;
end
