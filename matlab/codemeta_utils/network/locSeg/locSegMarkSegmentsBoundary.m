%Given a segmentation membership, mark members that constitute the border
%of each segment
%function members=locSegMarkSegmentsBoundary(members,E)
%For each line in members, nodes that belong to the border of the segment
% are marked with a "2" instead of a "1"
function members=locSegMarkSegmentsBoundary(members,E)
NNodes=size(members,2);
for iNode=1:NNodes
    idxPresent=members(:,iNode)>0;
    if sum(idxPresent)>1
        members(idxPresent,iNode)=2;
    end
end
        
% NMembers=size(members,1);
% for iMember=1:NMembers
%     for iNode=find(members(iMember,:))
%         NINode=findNeighbors(iNode,E);
%         if ~all(members(iMember,NINode)>0)
%             members(iMember,iNode)=2;
%         end
%     end
% end
        

