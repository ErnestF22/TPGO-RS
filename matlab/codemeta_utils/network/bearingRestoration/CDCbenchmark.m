addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codeMeta'))
addpath(genpath('C:\Users\Arman\Dropbox\GitHub\codeMeta'))
addpath(genpath('/Users/arman/Dropbox/GitHub/codeMeta'))
clc
clear
close all
sigmaNoise=0.0; tol=1e-12;
flagRemoveExcessive=true;

Niterations=5;

N=10:5:50;
Tave=cell(size(N));
Tvar=cell(size(N));

for Nidx=1:length(N)
    n=N(Nidx);
    IntervalI=n-1;
    IntervalF=floor(0.6*(n-1)*(n-2)/2);
    IntervalLength=IntervalF-IntervalI+1;
    TaveG=zeros(1,IntervalLength); TvarG=zeros(1,IntervalLength);

    for m=IntervalI:IntervalF
        disp([num2str(n),':',num2str(m)]);
        Tg=zeros(1,Niterations);
        Tc=zeros(1,Niterations);
        for iteration=1:Niterations
            %disp([num2str(m),':',num2str(iteration)]);
            resetRands(iteration-1);
            [E,x]=RandomGraphGenerator(n,m);
            u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
            membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
            [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
            EE=FactorGraph(Eprime,x,pins,membprime,nodemembership);

            tic
            [Erigid1,NaddedEdges,EE] = RigidifyGreedyBF(x,E);
            t = toc;
            Tg(iteration)=t;
            ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
            if ReqEdges~=NaddedEdges
                disp('DISCREPANCYgreedy');
            end       
        end
        TaveG(m-IntervalI+1)=mean(Tg);
        TvarG(m-IntervalI+1)=var(Tg);
        disp('---');
    end
    Tave(Nidx)={TaveG};
    Tvar(Nidx)={TvarG};
    hold on
    plot((IntervalI:IntervalF)*2/(n*(n-1)),TaveG,'DisplayName',['n=',num2str(n)]);
end

xlabel('D');
ylabel('Running Time (s)');
legend('show')
title(['#Iter=',num2str(Niterations)])