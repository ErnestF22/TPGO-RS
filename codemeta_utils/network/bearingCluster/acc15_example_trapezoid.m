function acc15_example_trapezoid
fileDir='../figures/tikz';
fileNameBase='trapezoid';


[A,x]=bearingCluster_generateTest('trapezoid');
E=adj2edges(A,'oriented');
u=bearingCluster_getBearingsScalesFromE(x,E);

membership=bearingCluster_clustering(E,u,'flagSeparateComponents',false);

x(1,[2 3])=x(1,[2 3])-0.5;

tikzOpts={'flagNormalizeScale',false,'membership',membership};
fileName=fullfile(fileDir,[fileNameBase '1.tex']);
fid=fopen(fileName,'wt');
bearingCluster_tikz(x,E,'fileId',fid,tikzOpts{:});
fclose(fid);

x(1,[2 3])=x(1,[2 3])+1.5;
fileName=fullfile(fileDir,[fileNameBase '2.tex']);
fid=fopen(fileName,'wt');
bearingCluster_tikz(x,E,'fileId',fid,tikzOpts{:});
fclose(fid);
