%function benchmark()
addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codeMeta'))
addpath(genpath('C:\Users\Arman\Dropbox\GitHub\codeMeta'))
addpath(genpath('/Users/armando/Dropbox/GitHub/codeMeta'))
clc
clear
close all
sigmaNoise=0.0; tol=1e-12;
flagRemoveExcessive=true;

n=7;%10
m=6;%11
resetRands(2);
[E,x]=RandomGraphGenerator(n,m);
%  E=[[1,2];[2,3];[3,4];[4,1]];
%  x=[[0;0],[1;0],[1;1],[0;1]];

MODE = 'Greedy'

switch MODE
    case 'Combinatorial'
        tic
        u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
        membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
        disp(['Number of Components: ',num2str(max(membership))]);
        [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
        ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
        [Erigid,NaddedEdges]=CombinatorialIndexing(E,x,membership,ReqEdges,'flagoperationmode','FirstOption');
        t = toc
        bearingClustering_plot(x,E,membership); figure;
        bearingClustering_plot(x,[E;Erigid]);
        disp('Combinatorial Done');
        if ReqEdges~=NaddedEdges
            disp('DISCREPANCYcomb');
        end
    case 'Greedy'
        tic
        u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
        membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
        [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
        EE=FactorGraph(Eprime,x,pins,membprime,nodemembership);
        %[RigidTable,NaddedEdges] = RigidifyGreedy(EE,E,membership,pins);
        [E,membership,NaddedEdges,EE] = RigidifyGreedy3d(EE,E,membership,pins,x);
        t = toc;
        ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
        if ReqEdges~=NaddedEdges
            disp('DISCREPANCYgreedy');
        end
        disp('Greedy Done');
end
 %bearingClustering_plot(x,E)