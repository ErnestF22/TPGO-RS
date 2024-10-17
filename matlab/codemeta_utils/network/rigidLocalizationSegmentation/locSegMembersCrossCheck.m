%Keeps only the rows that appear in both matrices
%function members=locSegMembersCrossCheck(members1,members2)
%NOTE: multiple copies of the same row in members1 will be kept, but
%multiple copies in members2 will be not kept. 

%TODO: maybe this function could be made more efficient
function [members,flagMember1]=locSegMembersCrossCheck(members1,members2)
NMembers1=size(members1,1);
NMembers2=size(members1,1);
flagMember1=zeros(1,NMembers1);
for iMember=1:NMembers1
    n=members1(iMember,:);
    %check if n is in members2
    if any(~any(members2-repmat(n,NMembers2,1),2))
        flagMember1(iMember)=1;
    end
end
members=members1(flagMember1,:);
