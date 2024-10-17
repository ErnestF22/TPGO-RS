%Merge memberships using internal/border labels for the nodes
%function members=locSegSplitJoinSegmentations3(members,newMembers)
function members=locSegSplitJoinSegmentations3(members,newMembers)
NNewMembers=size(newMembers,1);
NMembers=size(members,1);
membersAdd=[];
for iNewMember=1:NNewMembers
    nm=newMembers(iNewMember,:);
    for iMember=1:NMembers
        m=members(iMember,:);
        [mNew,mAdd]=splitJoinPair(m,nm);
        members(iMember,:)=mNew;
        membersAdd=[membersAdd;mAdd];
    end
end
members=[members;membersAdd];

function [mNew,mAdd]=splitJoinPair(m,nm)
%remove all elements from nm which are not in m
nm(m==0)=0;
if all(nm==2 | nm==0)
    %it is a compatible component
    mNew=m;
    mAdd=[];
else
    %split m
    %mIn is the part "inside" nm
    %mOut is the part "outside" nm
%     mIn=zeros(size(m));
%     mIn(nm==1)=1;
%     mIn(m==2 & nm==2)=2;
%     mOut=m;
%     mOut(nm==1)=0;
    mOut=zeros(size(m));
    mOut(m==1 & nm==0)=1;
    mOut(m==1 & nm==1)=0;
    mOut(m==1 & nm==2)=2;
    mOut(m==2 & nm==0)=2;
    mOut(m==2 & nm==1)=0;
    mOut(m==2 & nm==2)=2;
    mIn=zeros(size(m));
    mIn(m==1 & nm==0)=0;
    mIn(m==1 & nm==1)=1;
    mIn(m==1 & nm==2)=2;
    mIn(m==2 & nm==0)=0;
    mIn(m==2 & nm==1)=2;
    mIn(m==2 & nm==2)=2;
    %if all elements that are present in mOut are "2", just discard it
    if all(mOut==0 | mOut==2)
        mOut=[];
    end
    %if all elements in mIn where "2" originally, just discard it
    if all(m(mIn>0)==2) || all(mIn==0 | mIn==2)
        mIn=[];
    end
    %assign the non-empty element to mNew and, if present, the second one
    %to mAdd
    if ~isempty(mIn)
        mNew=mIn;
        mAdd=mOut;
    else
        mNew=mOut;
        mAdd=[];
    end
end
