function members=locSegAddSegmentToPool(members,newMembers)
NNewMembers=size(newMembers,1);
NMembers=size(members,1);

for iNewMember=1:NNewMembers
    n=newMembers(iNewMember,:);
    if isempty(members)
        %if members is empty, just add the row
        members=n;
    else
        %check wether n is already equal to one of the rows in members
        if ~any(~any(members-ones(NMembers,1)*n,2))
            %add row to members
            members=[members; n];
            NMembers=NMembers+1;
        end
    end
end
