function structOut=sfm_rawFilterDataWithFlag(structIn,flag,memberList)
if ~exist('memberList','var')
    memberList=fieldnames(structIn);
end

NDataExpected=length(flag);
NMember=length(memberList);
for iMember=1:NMember
    memberName=memberList{iMember};
    if size(structIn.(memberName),2)==NDataExpected
        structOut.(memberName)=structIn.(memberName)(:,flag);
    else
        structOut.(memberName)=structIn.(memberName);
    end
end
