%data set fot test_control function
function [tree,obstacles] = testData
tree(1).position = [25;28];
tree(2).position = [30;45];
tree(3).position = [35;20];
tree(4).position = [10;30];
tree(5).position = [20;10];
tree(6).position = [38;48];
tree(7).position = [25;50];
tree(8).position = [30;58];

tree(1).parent = 3;
tree(2).parent = 1;
tree(3).parent = [];
tree(4).parent = 1;
tree(5).parent = 1;
tree(6).parent = 2;
tree(7).parent = 2;
tree(8).parent = 7;

obstacles = [30 15 14 35;20 40 25 45];
end