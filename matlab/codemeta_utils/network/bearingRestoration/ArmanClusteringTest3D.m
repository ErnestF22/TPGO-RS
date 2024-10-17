%function ArmanClusteringTest3D
addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codemeta'))
addpath(genpath('../../'))
clc
clear
close all
%sigmaNoise=0.1; tol=0.1;
sigmaNoise=0.0; tol=1e-12;
%methodMeasurements='incidence';
methodMeasurements='cycleBasis';
flagSplitInterdependent=true;
flagRemoveExcessive=false;

c='RandomGraph';
 switch c
     case 'RandomGraph'
           resetRands(7)
           n=20;
           m=25;
         [E,x]=RandomGraphGenerator(n,m,'dimension',3);
    x=x+0.1*randn(size(x));
     case 'Triangle'
         E = [1,2;2,3;3,1];
         x = [0,0,0; 10,0,0; 10,10,0]';
     case 'RectanglePlanar'
         E = [1,2;2,3;3,4;1,4];
         x = [0,0,0; 10,0,0; 10,10,0; 0,10,0]';
     case 'Rectangle'
         E = [1,2;2,3;3,4;1,4];
         x = [0,0,0; 10,0,0; 10,10,0; 0,10,0.5]';
     case 'Pentagon'
         E = [1,2;2,3;3,4;4,5;1,5];
         x = [0,0,0.5; 10,0,0; 13,10,0; 10,13,0; 0,10,0]';
     case 'PentagonPlanar'
         E = [1,2;2,3;3,4;4,5;1,5];
         x = [0,0,0; 10,0,0; 13,10,0; 10,13,0; 0,10,0]';
     case 'PentagoalCycles'
         E=[1,2;2,3;3,4;4,5;1,5];
         E=[E; 1,6;6,7;7,8;2,8];
         E=[E; 8,9;9,10;3,10];
         E=[E; 10,11;11,12;4,12];
         E=[E; 12,13;13,14;5,14];
         E=[E; 14,15;6,15];
         E=[E; 1,7];
         x=zeros(3,15);
         L1=4;
         L2=7;
         for idx=1:5
             x(1,idx)=L1*cos(0.4*pi*idx);
             x(2,idx)=L1*sin(0.4*pi*idx);
             x(1,5+2*idx-1)=L2*cos(0.4*pi*idx);
             x(2,5+2*idx-1)=L2*sin(0.4*pi*idx);
             x(1,5+2*idx)=L2*cos(0.4*pi*idx+pi*0.2);
             x(2,5+2*idx)=L2*sin(0.4*pi*idx+pi*0.2);
         end
         x=x+0.2*randn(size(x));
     case 'Planar_Cycles'
         E=[1,2; 2,3; 3,4; 1,4];
         E=[E; 1,5; 2,5];
         x=zeros(3,max(E(:)));
         x(1:2,:)=randn(2,max(E(:)));
         x(3,5)=1;
 end

[u,lambda]=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);

[membership,L,s]=bearingCluster_clustering(E,u,'flagfirstcall',true);

%This code finds pins and removes excessive edges
[pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);

ReqEdges = ceil((3*size(x,2)-4)/2)-sum(ceil((3*sum(nodemembership,1)-4)/2));

%Code for Factor Graph
[EE,XX]=FactorGraph(Eprime,x,pins,membprime,nodemembership);

%Code for rigidifying the graph
[E,membership,NaddedEdges,EE] = RigidifyGreedy3d(EE,E,membership,pins,x);

fprintf('Framework contains:\n')
fprintf('\t%d shakes\n',size(L,2));
fprintf('\t%d rigid components\n',length(unique(membership)))
fprintf('\t%d fundamental cycles\n',size(EE,1)-size(XX,2)+1);
fprintf('\t%d pins\n',sum(pins));
fprintf('\t%d req edges\n',ReqEdges);
bearingClustering_plot(x,E,membership,'vertexNumber')
v = axis;


figure
FactorGraph_plot(XX,EE);
%axis(1.1*v);


% [Erigid,membership,NaddedEdges,EE] = RigidifyGreedy3d(EE,E,membership,pins,x)
% [u,lambda]=bearingCluster_getBearingsScalesFromE(x,Erigid,'noisy',sigmaNoise);
% [membership,L,s]=bearingCluster_clustering(Erigid,u,'flagfirstcall',true);
% membership
% 
% s1=getEigVal(Erigid,u);
% 
% 
% u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
% membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
% [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
% ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
% [Er,NaddedEdges]=CombinatorialIndexing(E,x,membership,ReqEdges);
% Erigid=[E;Er];
% [u,lambda]=bearingCluster_getBearingsScalesFromE(x,Erigid,'noisy',sigmaNoise);
% [membership,L,s]=bearingCluster_clustering(Erigid,u,'flagfirstcall',true);
% membership
% 
% s2=getEigVal(Erigid,u);
% 
% figure;
% hold on;
% plot(1:length(s1),s1,'-*r');
% plot(1:length(s2),s2,'-*k');
% legend('Greedy','Combinatorial');
% figure;
% bearingClustering_plot(x,Erigid,membership);
% %save([mfilename '_data'])
% %end