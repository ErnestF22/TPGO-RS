clc
clear
close all
addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codeMeta'))
addpath(genpath('C:\Users\Arman\Dropbox\GitHub\codeMeta'))
addpath(genpath('/Users/arman/Dropbox/GitHub/codeMeta'))
sigmaNoise=0.0; tol=1e-6;

x=zeros(2,8);
for nidx=1:8
    x(:,nidx)=[cos((2*nidx-1)*pi/8);sin((2*nidx-1)*pi/8)];
end

E=[1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1;1,4;5,8;2,7;3,6];
%E=[1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1;1,4;5,8;2,7;3,6;3,5];

u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
membership=bearingCluster_clustering(E,u,'flagfirstcall',true,'flagSeparateComponents',false);
membership2=bearingCluster_clustering(E,u,'flagfirstcall',true,'flagSeparateComponents',true);
[pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership2,'flagremoveexcessive',true);
ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
EE=FactorGraph(Eprime,x,pins,membprime,nodemembership);
bearingClustering_plot(x,E,membership,'vertexnumber');
grid off
axis(1.2*axis)

figure
bearingClustering_plot(x,E,membership2,'vertexnumber');
grid off
axis(1.2*axis)

[Erigid,NaddedEdges]=CombinatorialIndexing(E,x,membership2,ReqEdges,'flagoperationmode','BestL2');
figure
Er=[E;Erigid];
bearingClustering_plot(x,Er,ones(1,size(Er,1)),'vertexnumber');
grid off
axis(1.2*axis)