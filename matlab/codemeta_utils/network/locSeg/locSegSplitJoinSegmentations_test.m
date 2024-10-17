function locSegSplitJoinSegmentations_test
members=[
    1 1 0 0 0 0 0 0 0;
    0 0 1 1 0 0 0 0 0;
    0 0 0 0 0 1 1 1 1;
    ];
attributes={'a1','a2','a3'};

newMembers=[
    1 1 0 0 0 0 0 0 0;
    1 1 1 1 0 0 0 0 0;
    0 0 0 0 0 1 1 0 0;
];

newAttributes={'b1','b2','b3'};

disp(members)
disp(attributes)

disp('Adding same...')
[members,attributes]=locSegSplitJoinSegmentations(members,attributes,newMembers(1,:),newAttributes(1));

disp(members)
disp(attributes)

disp('Adding superset...')
[members,attributes]=locSegSplitJoinSegmentations(members,attributes,newMembers(2,:),newAttributes(2));

disp(members)
disp(attributes)

disp('Adding subset...')
[members,attributes]=locSegSplitJoinSegmentations(members,attributes,newMembers(3,:),newAttributes(3));

disp(members)
disp('The following should be {''a1'',''b1'',''b2''}')
disp(attributes{1})
disp('The following should be {''a2'',''b2''}')
disp(attributes{2})
disp('The following should be {''a3'',''b3''}')
disp(attributes{3})
disp('The following should be {''a3''}')
disp(attributes{4})


