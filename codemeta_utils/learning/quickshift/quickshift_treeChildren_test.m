%Test for the treeChildren function. 
%The expected output is:
% 1 ->  1	[2 6 ]
% 2 ->  1	[3 4 7 ]
% 3 ->  2	[5 ]
% 4 ->  2	[ ]
% 5 ->  3	[ ]
% 6 ->  1	[ ]
% 7 ->  2	[8 ]
% 8 ->  7	[ ]
% 9 ->  9	[10 ]
%10 ->  9	[11 ]
%11 -> 10	[ ]
%Note that every time node x is a children to y (i.e., x -> y), x appears
%in the list of children of y
function quickshift_treeChildren_test
treeEdges=[1;1;2;2;3;1;2;7;9;9;10];
treeChildren=quickshift_treeChildren(treeEdges);

for iEdge=1:length(treeEdges)
    fprintf('%2d -> %2d\t[', iEdge, treeEdges(iEdge))
    fprintf('%d ',treeChildren{iEdge})
    fprintf(']\n')
end
