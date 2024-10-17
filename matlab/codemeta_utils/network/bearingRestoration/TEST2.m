clc
clear
close all

sigmaNoise=0.0; tol=1e-12;
addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codeMeta'))
addpath(genpath('/Users/arman/Dropbox/GitHub/codeMeta'))
addpath(genpath('/Users/arman/Dropbox/GitHub/Code/CDC2017/SetPartFolder'))

n=8;
m=11;
iteration=8;
% n=10;
% m=14;
% iteration=21;
resetRands(iteration-1);
[E,x]=RandomGraphGenerator(n,m);

u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
membership=bearingCluster_clustering(E,u,'flagfirstcall',true,'threshold',0.2);
bearingClustering_plot(x,E,membership,'vertexnumber');
membership
figure
[Erigid,NaddedEdges,EE] = RigidifyGreedyBF(x,E);
u=bearingCluster_getBearingsScalesFromE(x,Erigid,'noisy',sigmaNoise);
membership=bearingCluster_clustering(Erigid,u,'flagfirstcall',true,'threshold',0.2);
bearingClustering_plot(x,Erigid,membership,'vertexnumber');

% [pins,nodemembership,Eprime,membprime] =
% PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
% EE=FactorGraph(Eprime,x,pins,membprime,nodemembership);
% [Erigid,NaddedEdges,EE] = RigidifyGreedyBF(x,E);
% 
% 
% 
% [R,S,flagIsRigid] = distanceRigidityMatrix(E,x);