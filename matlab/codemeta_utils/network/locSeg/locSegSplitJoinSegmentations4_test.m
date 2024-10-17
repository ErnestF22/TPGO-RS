function locSegSplitJoinSegmentations4_test
dataset=2;
switch dataset
    case 1
        E=[
             1     2     3     3     3     3     4     5     6     6     7     2     3     1     4     5     6     6     6     7     8     8
             2     3     1     4     5     6     6     6     7     8     8     1     2     3     3     3     3     4     5     6     6     7
        ]';
        membersEdge{1}=[
             0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1     1     1
             0     0     1     1     1     1     1     1     0     0     0     0     0     1     1     1     1     1     1     0     0     0
             1     1     1     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0
        ];
        membersEdge{2}=[
             0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1     1     1
             0     1     0     1     1     1     1     1     0     0     0     0     1     0     1     1     1     1     1     0     0     0
             1     1     1     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0
        ];
        membersEdge{3}=[
             0     0     0     0     0     1     0     0     1     1     1     0     0     0     0     0     1     0     0     1     1     1
             0     0     0     0     1     1     0     1     0     0     0     0     0     0     0     1     1     0     1     0     0     0
             0     0     0     1     0     1     1     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0
             1     1     1     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0
        ];
        membersEdge{4}=[
             0     0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0     1     0     1     1     1
             0     0     0     1     1     1     1     1     0     0     0     0     0     0     1     1     1     1     1     0     0     0
             1     1     1     1     0     0     0     0     0     0     0     1     1     1     1     0     0     0     0     0     0     0
        ];
        membersEdge{5}=[
             0     0     0     0     0     0     0     1     1     1     1     0     0     0     0     0     0     0     1     1     1     1
             0     0     0     1     1     1     1     1     0     0     0     0     0     0     1     1     1     1     1     0     0     0
             1     1     1     0     1     0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0     0
        ];
        membersEdge{6}=[
             0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1     1     1
             0     0     0     0     1     1     0     1     0     0     0     0     0     0     0     1     1     0     1     0     0     0
             0     0     0     1     0     1     1     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0
             1     1     1     0     0     1     0     0     0     0     0     1     1     1     0     0     1     0     0     0     0     0
        ];
        membersEdge{7}=[
             0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1     1     1
             0     0     0     1     1     1     1     1     1     0     0     0     0     0     1     1     1     1     1     1     0     0
             1     1     1     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0
        ];
        membersEdge{8}=[
             0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1     1     1
             0     0     0     1     1     1     1     1     0     1     0     0     0     0     1     1     1     1     1     0     1     0
             1     1     1     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0
        ];
    case 2
        E=[
     1     2     3     3     3     3     4     5     6     6     7     3     5     6     2     3     1     4     5     6     6     6     7     8     8     9     9     9
     2     3     1     4     5     6     6     6     7     8     8     9     9     9     1     2     3     3     3     3     4     5     6     6     7     3     5     6
        ]';
        membersEdge{1}=[
     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0
     0     0     1     1     1     1     1     1     0     0     0     1     1     1     0     0     1     1     1     1     1     1     0     0     0     1     1     1
     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{2}=[
     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0
     0     1     0     1     1     1     1     1     0     0     0     1     1     1     0     1     0     1     1     1     1     1     0     0     0     1     1     1
     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{3}=[
     0     0     0     0     0     1     0     0     1     1     1     0     0     0     0     0     0     0     0     1     0     0     1     1     1     0     0     0
     0     0     0     0     1     1     0     1     0     0     0     1     1     1     0     0     0     0     1     1     0     1     0     0     0     1     1     1
     0     0     0     1     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0     0     0     0
     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{4}=[
     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0     0
     0     0     0     1     1     1     1     1     0     0     0     1     1     1     0     0     0     1     1     1     1     1     0     0     0     1     1     1
     1     1     1     1     0     0     0     0     0     0     0     0     0     0     1     1     1     1     0     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{5}=[
     0     0     0     0     0     0     0     1     1     1     1     0     0     0     0     0     0     0     0     0     0     1     1     1     1     0     0     0
     0     0     0     1     1     1     1     1     0     0     0     1     1     1     0     0     0     1     1     1     1     1     0     0     0     1     1     1
     1     1     1     0     1     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{6}=[
     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0
     0     0     0     0     1     1     0     1     0     0     0     1     1     1     0     0     0     0     1     1     0     1     0     0     0     1     1     1
     0     0     0     1     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0     0     0     0
     1     1     1     0     0     1     0     0     0     0     0     0     0     0     1     1     1     0     0     1     0     0     0     0     0     0     0     0
        ];
        membersEdge{7}=[
     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0
     0     0     0     1     1     1     1     1     1     0     0     1     1     1     0     0     0     1     1     1     1     1     1     0     0     1     1     1
     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{8}=[
     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0
     0     0     0     1     1     1     1     1     0     1     0     1     1     1     0     0     0     1     1     1     1     1     0     1     0     1     1     1
     1     1     1     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0
        ];
        membersEdge{9}=[
     0     0     0     0     0     0     0     0     1     1     1     0     0     1     0     0     0     0     0     0     0     0     1     1     1     0     0     1
     0     0     0     1     1     1     1     1     0     0     0     1     1     1     0     0     0     1     1     1     1     1     0     0     0     1     1     1
     1     1     1     0     0     0     0     0     0     0     0     1     0     0     1     1     1     0     0     0     0     0     0     0     0     1     0     0
        ];
end        
stringsEdge=memberships2strings(membersEdge);
matches=stringsMatch(stringsEdge);
NEdges=size(E,1);
disp([zeros(2) E(1:NEdges/2,:)'; E(1:NEdges/2,:) matches(1:NEdges/2,1:NEdges/2)])
disp(matches)
C=bierstone(matches);
NC=length(C);
membersEdge=zeros(NC,NEdges);
for iC=1:NC
    membersEdge(iC,C{iC})=1;
end
disp([E';membersEdge])

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


function stringsEdge=memberships2strings(membersEdge)
NNodes=length(membersEdge);
NEdges=size(membersEdge{1},2);
stringsEdge=cell(1,NEdges);
for iEdge=1:NEdges
    stringsEdge{iEdge}=cellfun(@(x) find(x(:,iEdge)),membersEdge,'UniformOutput',false);
end
