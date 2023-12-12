%this script is a test function for checking the controller functions
close all
clear all

figure(2)
[tree,obs] = testData;
plotTree(tree,obs)
tree = ENVconstraint(tree);
tree = nearestObstacle(tree, obs);
tree = CLFconstraint(tree);
tree = CBFconstraint(tree, obs);
tree = findController(tree);
plotControllers(tree);



