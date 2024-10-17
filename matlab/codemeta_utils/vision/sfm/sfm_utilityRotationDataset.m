%Returns a dataset for absolute rotation estimation from relative measurements
%function [E,Rij,RiRef,RijNoise]=sfm_utilityRotationDataset()
function [E,Rij,RiRef,RijNoise]=sfm_utilityRotationDataset()
t_node=testNetworkBuildTestNetwork();
t_node=testNetworkAddMeasurements(t_node,'Method','Noisy');
t_node=splitgij(t_node);

E=t_node.E;
Rij=t_node.Rijtruth;
RiRef=t_node.Ritruth;
RijNoise=t_node.Rij;

[E,idxE]=edges2edges(E,'toOriented');
Rij=Rij(:,:,idxE);
RijNoise=RijNoise(:,:,idxE);