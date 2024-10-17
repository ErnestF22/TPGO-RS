addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codeMeta'))
addpath(genpath('C:\Users\Arman\Dropbox\GitHub\codeMeta'))
addpath(genpath('/Users/arman/Dropbox/GitHub/codeMeta'))
clc
clear
close all
sigmaNoise=0.0; tol=1e-12;
flagRemoveExcessive=false;
% 8:12:59 problem for BFS
% 12:18:41 problem for Brute Force Greedy
Niterations=5;

N=10:1:20;
Tave1=cell(size(N));
Tvar1=cell(size(N));
Tave2=cell(size(N));
Tvar2=cell(size(N));
Tave3=cell(size(N));
Tvar3=cell(size(N));
Lave1=cell(size(N));
Lvar1=cell(size(N));
Lave2=cell(size(N));
Lvar2=cell(size(N));
Lave3=cell(size(N));
Lvar3=cell(size(N));

for Nidx=1:length(N)
    n=N(Nidx);
    disp(n);
    IntervalI=n-1;
    IntervalF=floor(0.8*(n-1)*(n-2)/2);
    IntervalLength=IntervalF-IntervalI+1;
    TaveG=zeros(1,IntervalLength); TvarG=zeros(1,IntervalLength);
    TaveC1=zeros(1,IntervalLength); TvarC1=zeros(1,IntervalLength);
    TaveC2=zeros(1,IntervalLength); TvarC2=zeros(1,IntervalLength);
    LaveG=zeros(1,IntervalLength); LvarG=zeros(1,IntervalLength);
    LaveC1=zeros(1,IntervalLength); LvarC1=zeros(1,IntervalLength);
    LaveC2=zeros(1,IntervalLength); LvarC2=zeros(1,IntervalLength);
    for m=IntervalI:IntervalF
        %disp([num2str(n),':',num2str(m)]);
        Tg=zeros(1,Niterations);
        Tc1=zeros(1,Niterations);
        Tc2=zeros(1,Niterations);
        Lg=zeros(1,Niterations);
        L21=zeros(1,Niterations);
        L22=zeros(1,Niterations);
        for iteration=1:Niterations
            disp([num2str(n),':',num2str(m),':',num2str(iteration)]);
            resetRands(iteration-1);
            [E,x]=RandomGraphGenerator(n,m);
            u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
            membership=bearingCluster_clustering(E,u,'flagfirstcall',true,'threshold',0.2);
            [~,nodemembership,~,~] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
%             EE=FactorGraph(Eprime,x,pins,membprime,nodemembership);
            tic
            [Erigid1,NaddedEdges,EE] = RigidifyGreedyBF(x,E);
            t = toc;
            Tg(iteration)=t;
            u=bearingCluster_getBearingsScalesFromE(x,Erigid1,'noisy',sigmaNoise);
            eigval=getEigVal(Erigid1,u);
            Lg(iteration)=eigval(2);
            ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
            if ReqEdges~=NaddedEdges
                disp('DISCREPANCYgreedy');
%                 disp([num2str(n),':',num2str(m),':',num2str(iteration)]);
            end
            
            
            disp('---');
            
            tic
            u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
            membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
            [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
            ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
            [Erigid2,NaddedEdges,lambda2,Extra]=CombinatorialIndexing(E,x,membership,ReqEdges,'flagoperationmode','Both');
            t = toc;
            Tc1(iteration)=Extra(1);
            Tc2(iteration)=t;
            L21(iteration)=Extra(2);
            L22(iteration)=lambda2;
            if ReqEdges~=NaddedEdges
                disp('DISCREPANCYcomb');
            end
        end
        TaveG(m-IntervalI+1)=mean(Tg);
        TvarG(m-IntervalI+1)=var(Tg);
        TaveC1(m-IntervalI+1)=mean(Tc1);
        TvarC1(m-IntervalI+1)=var(Tc1);
        TaveC2(m-IntervalI+1)=mean(Tc2);
        TvarC2(m-IntervalI+1)=var(Tc2);
        LaveG(m-IntervalI+1)=mean(Lg);
        LvarG(m-IntervalI+1)=var(Lg);
        LaveC1(m-IntervalI+1)=mean(L21);
        LvarC1(m-IntervalI+1)=var(L21);
        LaveC2(m-IntervalI+1)=mean(L22);
        LvarC2(m-IntervalI+1)=var(L22);
    end
    
    
    Tave1(Nidx)={TaveG};
    Tvar1(Nidx)={TvarG};
    Tave2(Nidx)={TaveC1};
    Tvar2(Nidx)={TvarC1};
    Tave3(Nidx)={TaveC2};
    Tvar3(Nidx)={TvarC2};
    Lave1(Nidx)={LaveG};
    Lvar1(Nidx)={LvarG};
    Lave2(Nidx)={LaveC1};
    Lvar2(Nidx)={LvarC1};
    Lave3(Nidx)={LaveC2};
    Lvar3(Nidx)={LvarC2};
     
     figure(1)
     hold on
     plot((IntervalI:IntervalF)*2/(n*(n-1)),TaveG,'DisplayName',['n=',num2str(n)]);
%      plot((IntervalI:IntervalF)*2/(n*(n-1)),TaveG,'DisplayName',['n=',num2str(n), 'Greedy']);
%      plot((IntervalI:IntervalF)*2/(n*(n-1)),TaveC1,'DisplayName',['n=',num2str(n), 'Comb, First']);
%      plot((IntervalI:IntervalF)*2/(n*(n-1)),TaveC2,'DisplayName',['n=',num2str(n), 'Comb, Best']);
     
     figure(2)
     hold on
     plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveG,'DisplayName',['n=',num2str(n)]);
%      plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveG,'DisplayName',['n=',num2str(n), 'Greedy']);
%      plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveC1,'DisplayName',['n=',num2str(n), 'Comb, First']);
%      plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveC2,'DisplayName',['n=',num2str(n), 'Comb, Best']);
end


figure(1)
legend('show')
xlabel('D');
ylabel('Running Time (s)');
figure(2)
legend('show')
xlabel('D');
ylabel('\lambda_2');
% xlabel('D');
% ylabel('Running Time (s)');
% legend('show')
% title(['#Iter=',num2str(Niterations)])


% 
% %function benchmark()
% addpath(genpath('C:\Users\armandok\Dropbox\GitHub\codeMeta'))
% addpath(genpath('C:\Users\Arman\Dropbox\GitHub\codeMeta'))
% addpath(genpath('/Users/arman/Dropbox/GitHub/codeMeta'))
% clc
% clear
% close all
% sigmaNoise=0.0; tol=1e-12;
% flagRemoveExcessive=true;
% 
% Niterations=5;
% n=12;
% %m=120;
% IntervalI=11;
% IntervalF=50;
% IntervalLength=IntervalF-IntervalI+1;
% TaveG=zeros(1,IntervalLength); TvarG=zeros(1,IntervalLength);
% TaveC=zeros(1,IntervalLength); TvarC=zeros(1,IntervalLength);
% 
% for m=IntervalI:IntervalF
%     Tg=zeros(1,Niterations);
%     Tc2=zeros(1,Niterations);
%     for iteration=1:Niterations
%         disp([num2str(m),':',num2str(iteration)]);
%         resetRands(iteration-1);
%         [E,x]=RandomGraphGenerator(n,m);
%         tic
%         u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
%         membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
%         [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
%         EE=FactorGraph(Eprime,x,pins,membprime,nodemembership);
%         %[RigidTable,NaddedEdges] = RigidifyGreedy(EE,E,membership,pins);
%         [Erigid1,membership,NaddedEdges,EE] = RigidifyGreedy3(EE,E,membership,pins,x);
%         t = toc;
%         Tg(iteration)=t;
%         ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
%         if ReqEdges~=NaddedEdges
%             disp('DISCREPANCYgreedy');
%         end
%         disp('---');
%         tic
%         u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
%         membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
%         [pins,nodemembership,Eprime,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
%         ReqEdges = 2*size(x,2)-3-sum(2*sum(nodemembership)-3);
%         [Erigid2,NaddedEdges]=CombinatorialIndexing(E,x,membership,ReqEdges,'flagoperationmode','FirstOption');
%         t = toc;
%         Tc2(iteration)=t;
%         if ReqEdges~=NaddedEdges
%             disp('DISCREPANCYcomb');
%         end
%     end
%     TaveG(m-IntervalI+1)=mean(Tg);
%     TvarG(m-IntervalI+1)=var(Tg);
%     TaveC(m-IntervalI+1)=mean(Tc2);
%     TvarC(m-IntervalI+1)=var(Tc2);
% end
% 
% hold on
% plot(IntervalI:IntervalF,TaveG,'r')
% plot(IntervalI:IntervalF,TaveC,'k')
% legend('Greedy','Combinatorial')
% xlabel('m');
% ylabel('Running Time (s)');
% title(['n=',num2str(n),' #Iter=',num2str(Niterations)])