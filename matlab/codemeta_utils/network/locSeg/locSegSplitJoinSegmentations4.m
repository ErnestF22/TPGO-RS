function [membersFinal,C,matches]=locSegSplitJoinSegmentations4(members)
strings=memberships2strings(members);
matches=stringsMatch(strings);
C=bierstone(matches);

NC=length(C);
N=size(members{1},2);
membersFinal=zeros(NC,N);
for iC=1:NC
    membersFinal(iC,C{iC})=1;
end

function matches=stringsMatch(stringsEdge)
NEdges=length(stringsEdge);
matches=zeros(NEdges);
for iEdge=1:NEdges
    for jEdge=1:NEdges
        matches(iEdge,jEdge)=stringsMatchPair(stringsEdge{iEdge},stringsEdge{jEdge});
    end
end

function flag=stringsMatchPair(string1,string2)
l=length(string1);
flag=true;
for iL=1:l
    if isempty(intersect(string1{iL},string2{iL}))
        flag=false;
    end
end


function strings=memberships2strings(members)
N=max(cellfun(@(x) size(x,2),members));
idxNonEmptyMembers=find(~cellfun(@(x) isempty(x),members));
strings=cell(1,N);
for iN=1:N
    for iMember=idxNonEmptyMembers
        idx=find(members{iMember}(:,iN));
        if isempty(idx)
            idx=-1;
        end
        strings{iN}=[strings{iN} {idx}];
    end
end
