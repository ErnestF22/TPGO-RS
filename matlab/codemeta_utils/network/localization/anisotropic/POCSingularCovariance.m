function POCSingularCovariance
flagDisplay=false;

%resetRands(3)

methodAbsolutePoses='reference';

vBase1=[3;-3;0]+0.5*randn(3,1);

RBase=rot_randn(eye(3),0.1);
TBase=0.1*randn(3,1);
GBase=RT2G(RBase,TBase);

R1=rot_randn(eye(3),0.5);
T1=TBase+vBase1;
T1p=TBase+1.5*vBase1;
G1=RT2G(R1,T1);
G1p=RT2G(R1,T1p);

N=sphere_randn([0;0;1],0.5);
d=5;
NPoints=30;

NVec=planeNdToNVec(N,d);
X=planeGeneratePoints(NVec,NPoints);

[RPlane,TPlane]=planeToRT(NVec,'methodAbsolutePoses',methodAbsolutePoses);
NOrth=orthComplement(N);
vPlanep=3*NOrth*randn(2,1);
TPlanep=TPlane+vPlanep;
RPlanep=RPlane*rot(rand*[0;0;1]);
GPlane=RT2G(RPlane,TPlane);
GPlanep=RT2G(RPlanep,TPlanep);

GBase1=computeRelativePoseFromG(GBase,G1,'methodAbsolutePoses',methodAbsolutePoses);
GBase1p=computeRelativePoseFromG(GBase,G1p,'methodAbsolutePoses',methodAbsolutePoses);
SigmaBase1=computeSigmaEpipole(GBase1);
disp(computeError(GBase1,GBase1p,SigmaBase1))

GBasePlane=computeRelativePoseFromG(GBase,GPlane,'methodAbsolutePoses',methodAbsolutePoses);
GBasePlanep=computeRelativePoseFromG(GBase,GPlanep,'methodAbsolutePoses',methodAbsolutePoses);
SigmaBasePlane=computeSigmaPlane(GBasePlane);
disp(computeError(GBasePlane,GBasePlanep,SigmaBasePlane))


if flagDisplay
    figure(1)
    draw3dcameraFromG(GBase,'methodAbsolutePoses',methodAbsolutePoses)
    hold on
    draw3dcameraFromG(G1,'methodAbsolutePoses',methodAbsolutePoses)
    draw3dcameraFromG(G1p,'methodAbsolutePoses',methodAbsolutePoses,'color1','r','color2','r')
    draw3dcameraFromG(GPlane,'methodAbsolutePoses',methodAbsolutePoses)
    draw3dcameraFromG(GPlanep,'methodAbsolutePoses',methodAbsolutePoses,'color1','r','color2','r')
    plotPoints(X)
    hold off
    axis equal
end

function Sigma=computeSigmaEpipole(G)
T=G2T(G);
Sigma=blkdiag(eye(3),orthComplementProjector(T));

function Sigma=computeSigmaPlane(G)
R=G2R(G);
e3=[0;0;1];
Sigma=blkdiag(orthComplementProjector(e3),R(:,3)*R(:,3)');


function e=computeError(G,GTruth,Sigma)
v=rot3r3_vee(GTruth,rot3r3_log(GTruth,G));
e=(v'*Sigma*v)/2;
