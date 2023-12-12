function POCbreakTreeMergeDebug
load POCbreakTreeMergeDebug_dataInput.mat
quickshift_breakTreeMerge(treeEdges,treeData,edgesSorted,fMerge,'debug');
disp('Function returned without finding any inconsistency')
