function [C,Nodes]=BFS_cycle(EE,n0)
% clc
% clear
% EE=[1,6;3,6;3,5;2,5;2,4;1,4;1,7];
% n0=7;

global Nodemap NPins NNodes
NPins=max(EE(:,1));
Nodemap=sort(unique(EE(:)));
NNodes=length(Nodemap);
NEdges=size(EE,1);

if NEdges-NNodes+1==0   %IF there is no cycle
    C=[];
    Nodes=[];
    return
end

Current=false(NNodes,1);
Next=false(NNodes,1);
Scanned=false(NNodes,1);
Parent=zeros(NNodes,1);

N0=Adr(n0);
Current(N0)=true;

D=0;
while 1
    K=find(Current,1); %Next node in the current queue (stage)
    if isempty(K)   %Current stage done
        D=D+1;  %Go to the next stage
        Current=Next;
        Next=false(NNodes,1);
        K=find(Current,1);
        if isempty(K)
            disp('No Cycles!');
            break
        end
    end
    Current(K)=0;
    List=Neighbors(Nodemap(K),EE);
    List=Adr(List);
    temp=false(NNodes,1);
    temp(List)=true;
    List=temp;
    if K~=N0
        List(Parent(K))=false;
    end
    Next=Next|List;
    Scanned=Scanned+List;
    CircPoint=find(Scanned==2,1);
    if ~isempty(CircPoint)
        %disp('Cycle Found!');
        Knew=find(Scanned==2,1); %WHAT IF MORE THAN ONE CYCLE DISCOVERED?
        Nodes=Knew;
        while Parent(Nodes(1))~=0
            Nodes=[Parent(Nodes(1)),Nodes];
        end
        Nodes=[Nodes,K];
        while Parent(Nodes(end))~=0
            Nodes=[Nodes,Parent(Nodes(end))];
        end
        index=1;
        if sum(Nodes==Nodes(index))==2
            while sum(Nodes==Nodes(index))==2
                index=index+1;
            end
            index=index-1;
        end
        bounds=find(Nodes==Nodes(index));
        Nodes=Nodes(bounds(1):bounds(2)-1);
        Nodes=Nodemap(Nodes);
        break
    else
        Parent(List)=K;
    end
end
L=numel(Nodes)/2;
C=zeros(2*L,1);
%C=zeros(1,NEdges);
for idx=1:2*L-1
    ee=sort([Nodes(idx),Nodes(idx+1)]);
    [~,indx]=ismember(EE,ee,'rows');
    C(idx)=find(indx,1);
    %C=C+indx.';
end
ee=sort([Nodes(1),Nodes(end)]);
[~,indx]=ismember(EE,ee,'rows');
C(idx+1)=find(indx,1);
%C=C+indx.';
end

function [val]=Adr(n) % Address of node n Nodemap -> 1:NNodes
global Nodemap
val=zeros(length(n),1);
for idx=1:length(n)
    val(idx)=find(Nodemap==n(idx));
end
end

function [val]=Neighbors(n,EE) %Neighbors of n
global NPins
%val=false(NNodes,1);
if n>NPins
    list=EE(EE(:,2)==n,1);
else
    list=EE(EE(:,1)==n,2);
end
%val(list)=true;
val=list;
end