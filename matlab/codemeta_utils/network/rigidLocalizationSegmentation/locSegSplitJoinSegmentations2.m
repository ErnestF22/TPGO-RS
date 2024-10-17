%Merge memberships using internal/border labels for the nodes
%function members=locSegSplitJoinSegmentations2(members,newMembers)
%If newMembers is empty, exit without touching members
%Otherwise, for each row in newMembers
%   - If is already present, exit without touching members
%   - If is compatible with all the existing members, add it and pass to
%     the next row 
%   - Otherwise, find incompatible rows, split them, and do a recursive
%     call to add them
function members=locSegSplitJoinSegmentations2(members,newMembers,E)
NNewMembers=size(newMembers);

for iNewMember=1:NNewMembers
    nm=newMembers(iNewMember,:);
    if ~isPresent(members,nm)
        if isCompatible(members,nm)
            members=[members;nm];
        else
            [members,membersAdd,nm]=splitJoin(members,nm,E);
            members=locSegSplitJoinSegmentations2(members,nm,E);
            members=locSegSplitJoinSegmentations2(members,membersAdd,E);
        end
    end
end

function [members,membersAdd,nm]=splitJoin(members,nm,E)
NMembers=size(members,1);
flagMembers=ones(size(members,1),1);
membersAdd=[];
for iMember=1:NMembers
    m=members(iMember,:);
    if all(nm==m)
        nm=[];
        break
    else
        if isCompatiblePair(nm,m)
            continue
        else
            %is not compatible and is not a repetition
            flagMembers(iMember)=0;
            membersAdd=[membersAdd; AminusB(m,nm,E); AintersectionB(m,nm)];
            nm=AminusB(nm,m,E);
        end
    end
    if isempty(nm)
        break
    end
end
members=members(flagMembers==1,:);


function flag=isPresent(members,m)
flag=false;
for iMembers=1:size(members,1)
    if all(members(iMembers,:)==m)
        flag=true;
        break
    end
end

function flag=isCompatible(members,m)
flag=true;
for iMembers=1:size(members,1)
    if ~isCompatiblePair(members(iMembers,:),m)
        flag=false;
        break
    end
end

function flag=isCompatiblePair(mA,mB)
flagOverlap=(mA.*mB)>0;
flag=all(flagOverlap==0) || (all(mA(flagOverlap)==2) && all(mB(flagOverlap)==2));

function m=AminusB(mA,mB,E)
m=zeros(size(mA));
flagA1=mA>0 & mB==0;
flagA2=mA==1 & mB==2;
flagA3=mA==2 & mB==2;
m(flagA1)=mA(flagA1);
m(flagA2)=2;
m(flagA3)=2;
m=pruneNonBorder(m,flagA2 | flagA3,E);
if all(m==0) || (all(flagA1==0) && all(flagA2==0))
    m=[];
end

function [mNew,flag]=pruneNonBorder(m,flag,E)
mNew=m;
flagNew=flag;
for iNode=find(flag)
    NINode=findNeighbors(iNode,E);
    %if all neighbors have been modified or are zero, they are either zero or two, with
    %at least one equal to two, then remove the node
    NNINode0=sum(m(NINode)==0);
    NNINode2=sum(m(NINode)==2);
    if all(flag(NINode)==1 | m(NINode)==0) && NNINode0+NNINode2==length(NINode) && NNINode2>0
        mNew(iNode)=0;
        flagNew(iNode)=0;
    end
end

function m=AintersectionB(mA,mB)
m=zeros(size(mA));
flagA1=mA==1 & mB==1;
flagA2=(mA==1 & mB==2) | (mA==2 & mB==1);
flagA3=mA==2 & mB==2;
m(flagA1)=1;
m(flagA2)=2;
m(flagA3)=2;
if all(m==0) || (all(flagA1==0) && all(flagA2==0))
    m=[];
end
