function rotLocConsensusTestTrials
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

funs=rot3_almostGlobal_functions('type','tron','b',3);

N=8;
E=adj2edges(adjgallery(N,'banded',2));
maxIt=5000;
epsilon=1/(2*edges2maxDegree(E)*funs.mumax);

NConfigs=5;
NTrials=10;
errors=cell(NConfigs,NTrials);
data=cell(NConfigs);

for iConfig=1:NConfigs
    fprintf('# Config %d/%d\n',iConfig,NConfigs)
    R=rot_randn(eye(3),[],N);
    RRel=computeRelativeRotations(E,R);
    data{iConfig}.R=R;
    data{iConfig}.RRel=RRel;
    for iTrial=1:NTrials
        fprintf('## Trial %d/%d\n',iTrial,NTrials)
        e.RInit=cat(3,eye(3),rot_randn(eye(3),1000,N-1));
        
        disp('### Frobenious')
        fprintf('Simulation: ')
        [RFrob,QFrob]=rotLocFrobConsensus(E,RRel,'collect','showProgress',...
            'maxIt',maxIt,'RInit',e.RInit);
        e.R.Frobenious=RFrob;
        e.Q.Frobenious=QFrob;
        fprintf('Cost:       ')
        e.c.Frobenious=rotLocFrobCost(E,QFrob,RRel,'showProgress');
        fprintf('Distances:  ')
        e.d.Frobenious=alignErrors(R,RFrob);
        disp(mean(e.d.Frobenious(end,:)))
        
        disp('### Riemannian')
        fprintf('Simulation: ')
        RRiemannian=rotLocRiemannianConsensus(E,RRel,'collect','showProgress',...
            'maxIt',maxIt,'funs',funs,'epsilon',epsilon,'RInit',e.RInit);
        e.R.Riemannian=RRiemannian;
        fprintf('Cost:       ')
        e.c.Riemannian=rotLocRiemannianCost(E,RRiemannian,RRel,funs,'showProgress');
        fprintf('Distances:  ')
        e.d.Riemannian=alignErrors(R,RRiemannian);
        disp(mean(e.d.Riemannian(end,:)))

        errors{iConfig,iTrial}=e;
        save(fileNameSave,'errors','data')
    end
end

save(fileNameSave)

function RRel=computeRelativeRotations(E,R)
NEdges=size(E,1);
RRel=zeros([size(R,1) size(R,2) NEdges]);
for iEdge=1:NEdges
    Ri=R(:,:,E(iEdge,1));
    Rj=R(:,:,E(iEdge,2));
    RRel(:,:,iEdge)=Ri'*Rj;
end

function e=alignErrors(R,REst)
NR=size(R,3);
NIt=size(REst,4);
e=zeros(NIt,NR);
w=getTextWaitBar(NIt);
w(0)
for it=1:NIt
    e(it,:)=rotationProcrustesAlignError(R,REst(:,:,:,it));
    w(it)
end
