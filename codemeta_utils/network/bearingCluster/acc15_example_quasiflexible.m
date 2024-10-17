function acc15_example_quasiflexible
resetRands(1)
s1=writeFile('3-connected',0);
s2=writeFile('3-connected-skew',0);
s3=writeFile('3-connected',0.2,0.3);

setFigFont('Times')
fs=setFigFontSize(8);
plot([s1 s2 s3],'.-','MarkerSize',10)
legend('Flexible','Perturbed','Noisy','Location','SouthWest')
axis([0 16 0 1.1])
savefigure('../figures/3-connected-singularvalues','epsc',[200 100],2)
setFigFontSize(fs)

writeFile('3-connected-skew',0,0.3);
writeFile('3-connected',0.2,0.3);

function s=writeFile(datasetName,sigmaNoise,tol)
if ~exist('sigmaNoise','var')
    sigmaNoise=0;
end
if ~exist('tol','var')
    tol=1e-12;
end
[A,x]=bearingCluster_generateTest(datasetName);
E=adj2edges(A,'oriented');
u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
[membership,L,s]=bearingCluster_clustering(E,u,'tol',tol);

fileDir='../figures/tikz';
fileNameBase=datasetName;
if sigmaNoise>0
    fileNameBase=[fileNameBase '_noisy'];
end
if tol>1e-10
    fileNameBase=[fileNameBase '_thresholded'];
end

fileName=fullfile(fileDir,[fileNameBase '.tex']);
bearingCluster_tikz(x,E,'fileId',fileName,'membership',membership,...
    'flagNormalizeScale',true,'scale',2.5);
