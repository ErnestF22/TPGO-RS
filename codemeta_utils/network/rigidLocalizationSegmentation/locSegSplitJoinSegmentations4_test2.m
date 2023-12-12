function locSegSplitJoinSegmentations4_test2
dataset=1;
switch dataset
    case 1
        S=[
     1     2     3     4     5     6     7     8     9
     ];
        members{1}=[
     1     0     0     0     0     1     1     1     0
     1     0     1     1     1     1     0     0     1
     1     1     1     0     0     0     0     0     0
        ];
        members{2}=[
     0     1     0     0     0     1     1     1     0
     0     1     1     1     1     1     0     0     1
     1     1     1     0     0     0     0     0     0
        ];
        members{3}=[
     0     0     1     0     0     1     1     1     0
     0     0     1     0     1     1     0     0     1
     0     0     1     1     0     1     0     0     0
     1     1     1     0     0     0     0     0     0
        ];
        members{4}=[
     0     0     0     1     0     1     1     1     0
     0     0     1     1     1     1     0     0     1
     1     1     1     1     0     0     0     0     0
        ];
        members{5}=[
     0     0     0     0     1     1     1     1     0
     0     0     1     1     1     1     0     0     1
     1     1     1     0     1     0     0     0     0
        ];
        members{6}=[
     0     0     0     0     0     1     1     1     0
     0     0     1     0     1     1     0     0     1
     0     0     1     1     0     1     0     0     0
     1     1     1     0     0     1     0     0     0
        ];
        members{7}=[
     0     0     0     0     0     1     1     1     0
     0     0     1     1     1     1     1     0     1
     1     1     1     0     0     0     1     0     0
        ];
        members{8}=[
     0     0     0     0     0     1     1     1     0
     0     0     1     1     1     1     0     1     1
     1     1     1     0     0     0     0     1     0
        ];
        members{9}=[
     0     0     0     0     0     1     1     1     1
     0     0     1     1     1     1     0     0     1
     1     1     1     0     0     0     0     0     1
        ];
end

[~,C,matches]=locSegSplitJoinSegmentations4(members);

fprintf('Cross-segmentation matches\n')
disp(matches)

fprintf('Combined segmentation\n')
cellfun(@(x) disp(x),C,'UniformOutput',false)

% NC=length(C);
% N=size(members{1},2);
% members2=zeros(NC,N);
% for iC=1:NC
%     members2(iC,C{iC})=1;
% end
% fprintf('Combined segmentation\n')
% disp([S;members2])

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
N=size(members{1},2);
strings=cell(1,N);
for iN=1:N
    strings{iN}=cellfun(@(x) find(x(:,iN)),members,'UniformOutput',false);
end
