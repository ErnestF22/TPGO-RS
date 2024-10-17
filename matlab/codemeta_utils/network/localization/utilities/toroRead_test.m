function toroRead_test
fileName='../../../datasets/poseGraph/INTEL_P_toro.graph';
t_node=toroRead(fileName);
testNetworkDisplayErrors(t_node,'rt')
