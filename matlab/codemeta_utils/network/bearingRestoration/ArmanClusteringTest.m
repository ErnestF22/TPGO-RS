%function ArmanClusteringTest
addpath(genpath('../../'))
%addpath(genpath('/Users/arman/Dropbox/GitHub/codemeta'))
%addpath(genpath('/Users/arman/Dropbox/codeMeta'))
clc
clear
close all
%sigmaNoise=0.1; tol=0.1;
sigmaNoise=0.0; tol=1e-12;
%methodMeasurements='incidence';
methodMeasurements='cycleBasis';
flagSplitInterdependent=true;
flagRemoveExcessive=true;

c='Cycle_around_pentagon';
 switch c
     case 'file'
         load('testData/testData1.mat','x','E')
     case 'preresults'
         load('data.mat','x','E')
     case 'DoubleCycle'
         x=[];
         E=[];
         for idx=1:5
             x=[x,[cosd(idx*72+36);sind(idx*72+36)],2*[cosd(idx*72+36);sind(idx*72+36)]];
             E=[E;[2*idx-1,2*idx]];
         end
         x=x+0.1*randn(size(x));
         E=[E;1,3;3,5;5,7;7,9;1,9;2,4;4,6;6,8;8,10;2,10];
         E=[E;1,4;3,6];
         %E=[E;1,5;3,7];
     case 'Rigid7'
         x=randn(2,7);
         E=[1,2;2,3;1,3;2,4;1,4;2,5;5,6;1,6;1,5;1,7;2,7];
     case 'RandomGraph'
           resetRands(7)
           n=10;
           m=20;
%            resetRands(1)
%            n=20;
%            m=21; 
%           resetRands(0)
%           n=5;
%           m=5;
         [E,x]=RandomGraphGenerator(n,m);    
     case 'butterfly'
         E = [1,2; 2,3; 3,4; 4,1; 2,4; 3,5; 5,6; 6,7; 7,3; 5,7; 1,6];
         x = [50,50; 100,50; 100,100; 50,100; 150,100; 150,150; 100,150]';
     case 'butterflyx'
         E = [1,2; 2,3; 3,4; 4,1; 2,4; 3,5; 5,6; 6,7; 7,3; 5,7; 1,6; 7,8];
         x = [50,50; 100,50; 100,100; 50,100; 150,100; 150,150; 100,150; 75,175]';
     case 'hourglass'
         E = [1,2; 2,3; 1,3; 3,4; 4,5; 3,5; 5,6; 6,2; 2,5; 1,8; 1,7; 7,8; 8,9; 9,10; 8,10; 10,7; 10,11; 7,11; 11,6; 4,9];
         x = [0,0; -5,5; 5,5; 7,5; 0,10; -7,5; 5,-5; -5,-5; -7,-5; 0,-10; 7,-5]';
     case 'disturbedhourglass'
         E = [1,2; 2,3; 1,3; 3,4; 4,5; 3,5; 5,6; 6,2; 2,5; 1,8; 1,7; 7,8; 8,9; 9,10; 8,10; 10,7; 10,11; 7,11; 11,6; 4,9];
         x = [0,0; -5,5; 5,5; 7,5; 0,10; -7,5; 2.5,-2.5; -2.5,-2.5; -3.5,-2.5; 0,-5; 3.5,-2.5]';
     case 'test'
         E = [2,3; 3,4; 4,5; 5,2; 3,5; 4,6; 6,7; 7,8; 8,4; 6,8; 2,7];
         x = [50,50; 100,50; 100,100; 50,100; 150,100; 150,150; 100,150]';
     case 'CrossedRectangles'
         E = [1,2;2,3;3,4;4,1;1,6;3,6;3,5;1,5];
         x = [0,10;10,10;10,0;0,0;-8,-8;18,18]';
         x=x+randn(size(x));
     case 'CrossedRectangles2'
         E = [1,2;2,3;3,4;4,1;1,6;3,6;3,5;1,5;1,7;7,3;3,8;8,1];
         x = [0,10;10,10;10,0;0,0;-8,-8;18,18;16,16;-6,-6]';
         x=x+randn(size(x));
     case 'PerfectOctagon'
         E = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1;1,6;2,5];
         a = sqrt(2)/2;
         x = [1,0;a,a;0,1;-a,a;-1,0;-a,-a;0,-1;a,-a]';
     case 'Rectangle'
         E = [1,2;2,3;3,4;4,1];
         x = [0,0; 10,0; 10,10; 0,10]';
     case 'Triangle'
         E = [1,2;2,3;3,1];
         x = [0,0; 10,0; 10,10]';
     case 'Petersen'
         R1=2;
         R2=4;
         E = [1,3; 1,4; 2,4; 2,5; 3,5];
         E = [E;   6,7; 7,8; 8,9; 9,10; 6,10];
         E = [E;   1,6; 2,7; 3,8; 4,9; 5,10];
         %E = [E; 1,2; 1,5];  %Comment this line if need actual Petersen
         E = [E; 1,2; 3,4];  %Comment this line if need actual Petersen
         x=zeros(2,10);
         for idx=1:5
            v1=R1*[cos(0.4*pi*idx);sin(0.4*pi*idx)];
            v2=R2*[cos(0.4*pi*idx);sin(0.4*pi*idx)];
            x(:,idx)=v1;
            x(:,5+idx)=v2;
         end
         x=x+0.2*randn(size(x));
     case 'Utility'
         E = [1,4; 1,5; 1,6; 2,4; 2,5; 2,6; 3,4; 3,5; 3,6];
         x = [[0;0], [2;0], [4;0], [0;3], [2;3], [4;3]];
         x=x+0.1*randn(size(x));
         %x=2*randn(size(x));
     case 'Cycle_around_pentagon'
         E=zeros(15,2);
         x=zeros(2,10);
         v=[1:5,1];
         r1=2; r2=3.5;
         for idx=1:5
             E(idx,:)=[v(idx),v(idx+1)];
             E(5+idx,:)=[5+v(idx),5+v(idx+1)];
             E(10+idx,:)=[v(idx),5+v(idx)];
             x(:,idx)=r1*[cos(0.4*pi*idx);sin(0.4*pi*idx)];
             x(:,idx+5)=r2*[cos(0.4*pi*idx);sin(0.4*pi*idx)];
         end
         %E=[E; 1,3; 3,5];
         E=[E; 5,6; 2,8];
         x=x+0.1*randn(size(x));
 end
%E = uint8(E);

[u,lambda]=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);

[membership,L,s]=bearingCluster_clustering(E,u,'flagfirstcall',true);

%This code finds pins and removes excessive edges
[pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);

ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership,1)-3);

%Code for Factor Graph
[EE,XX]=FactorGraph(Eprime,x,pins,membprime,nodemembership);

%Code for rigidifying the graph
%[RigidTable,NaddedEdges] = RigidifyGreedy(EE,membership,pins)

fprintf('Framework contains:\n')
fprintf('\t%d shakes\n',size(L,2));
fprintf('\t%d rigid components\n',length(unique(membership)))
fprintf('\t%d fundamental cycles\n',size(EE,1)-size(XX,2)+1);
fprintf('\t%d pins\n',sum(pins));
fprintf('\t%d req edges\n',ReqEdges);
bearingClustering_plot(x,E,membership,'vertexNumber')
v = axis;

if membership~=membprime
    figure;
    bearingClustering_plot(x,Eprime,membprime,'vertexNumber')
    axis(v);
end

figure
FactorGraph_plot(XX,EE);
%axis(1.1*v);


% [Erigid,membership,NaddedEdges,EE] = RigidifyGreedy(EE,E,membership,pins,x)
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