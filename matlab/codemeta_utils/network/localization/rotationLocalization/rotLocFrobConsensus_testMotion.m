function rotLocFrobConsensus_testMotion
resetRands()

%Graph and initial conditions
t_node=testNetworkBuildTestNetwork('references');
RTruth0=t_node.Ritruth;
E=t_node.E;
NNodes=size(RTruth0,3);
NEdges=size(E,1);
RInit=rot_randn(RTruth0,10);

%Trajectories of the rotations
T=12;
Nt=T*5;
t=linspace(0,T,Nt);
v0=0.05*[zeros(3,1) randn(3,NNodes-1)];
vt=@(t) thetaFun(t,T)*v0;
RTrutht=@(t) rot_expVec(RTruth0,vt(t));

disp('Mean angular velocity [deg]')
disp(rad2deg(mean(cnorm(v0(:,2:end)))))


%Relative poses
RTruth=zeros([3 3 NNodes Nt]);
RRel=zeros([3 3 NEdges Nt]);
for it=1:Nt
    RTruth(:,:,:,it)=RTrutht(t(it));
    RRel(:,:,:,it)=relativeRotations(RTruth(:,:,:,it),E,'edges');
end

%Estimation
NIter=[5, 10, 30];
NNIter=length(NIter);
REst=zeros(size(RTruth));
eAligned=zeros(NNodes,Nt,NNIter);
for iNIter=1:NNIter
    NIt=NIter(iNIter);
    fprintf('NIter: %d: ',NIt);
    w=getTextWaitBar(Nt);
    w(0)
    for iT=1:Nt
        if iT==1
            REst(:,:,:,iT)=rotLocFrobConsensus(E,RRel(:,:,:,iT),'maxIt',NIt);
        else
            REst(:,:,:,iT)=rotLocFrobConsensus(E,RRel(:,:,:,iT),'maxIt',NIt,...
                'RInit',REst(:,:,:,iT-1));
        end

        eAligned(:,iT,iNIter)=rotationProcrustesAlignError(REst(:,:,:,iT),RTruth(:,:,:,iT));
        w(iT)
    end
end

legendText=[repmat('# Iterations = ',NNIter,1) num2str(NIter','%2d')];
plot(t,rad2deg(squeeze(mean(eAligned)))')
legend(legendText)
xlabel('s')
ylabel('Mean error [deg]')

save([mfilename '_data'])


function a=thetaFun(t,T)
tBreak=T/4;
a=zeros(size(t));
idxFirst=t<tBreak;
a(idxFirst)=t(idxFirst);
tHold=tBreak;
idxSecond=and(t>=tBreak,t<2*tBreak);
a(idxSecond)=tHold;
tNext=2*tBreak;
idxThird=t>=2*tBreak;
a(idxThird)=t(idxThird)-tNext+tHold;
