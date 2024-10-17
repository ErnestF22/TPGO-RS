function locSegSplitJoinSegmentations2_test
dataset=3;
switch dataset
    case 1
        members1=[
     2     0     0     2     0     1 
     2     0     2     2     1     0 
     2     1     2     0     0     0 
        ];
        members2=[
     2     0     0     2     0     1
     0     0     2     2     1     0
     2     1     2     2     0     0
        ];
        members3=[
     2     0     2     2     0     1
     0     0     2     2     1     0
     2     1     2     0     0     0
        ];
        membersGoal=[
     2     1     2     0     0     0
     2     0     0     2     0     1
     2     0     2     2     0     0
     0     0     2     2     1     0
        ];
    case 2
        members1=[
     2     0     0     2     0     1     0     0
     2     0     2     2     1     0     0     1
     2     1     2     0     0     0     1     0
     ];

        members2=[
     2     0     0     2     0     1     0     0
     0     0     2     2     1     0     0     1
     2     1     2     2     0     0     1     0
     ];

        members3=[
     2     0     2     2     0     1     0     0
     0     0     2     2     1     0     0     1
     2     1     2     0     0     0     1     0
     ];

        membersGoal=[
     2     1     2     0     0     0     1     0
     2     0     0     2     0     1     0     0
     2     0     2     2     0     0     0     0
     0     0     2     2     1     0     0     1
     ];
    case 3
        E=[
            1,2;
            2,3;
            3,4;
            4,1
            ];
        members1=[
     2     0     2     1
     2     1     2     0
     ];
        members2=[
     1     2     0     2
     0     2     1     2
     ];
        members3=[
     2     0     2     1
     2     1     2     0
     ];
        membersGoal=[
     2     2     0     0
     0     2     2     0
     0     0     2     2
     2     0     0     2
     ];
            
 
 
        
end

fprintf('Sequence 1 2 3\n')
members=locSegSplitJoinSegmentations2(members1,members2,E);
members=locSegSplitJoinSegmentations2(members,members3,E);
checkResult(members,membersGoal)

fprintf('Sequence 1 3 2\n')
members=locSegSplitJoinSegmentations2(members1,members3,E);
members=locSegSplitJoinSegmentations2(members,members2,E);
checkResult(members,membersGoal)

fprintf('Sequence 3 2 1\n')
members=locSegSplitJoinSegmentations2(members3,members2,E);
members=locSegSplitJoinSegmentations2(members,members1,E);
checkResult(members,membersGoal)


function checkResult(members,membersGoal)
members=sortrows(members);
membersGoal=sortrows(membersGoal);
fprintf('members\n')
disp(members)
fprintf('membersGoal\n')
disp(membersGoal)
if numel(members)==numel(membersGoal)
    fprintf('Difference: %s\n', num2str(max(abs(members(:)-membersGoal(:)))))
else
    fprintf('Different number of elements\n')
end

