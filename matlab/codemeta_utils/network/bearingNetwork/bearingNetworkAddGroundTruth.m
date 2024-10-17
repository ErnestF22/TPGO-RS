%Add ground truth poses and measurements (bearing and, possibly, ranges)
%function t_node=bearingNetworkAddGroundTruth(t_node,Titruth)
%If Titruth is not provided, evenly distribute the nodes around a circle
function t_node=bearingNetworkAddGroundTruth(t_node,Titruth)
NNodes=t_node.NNodes;

if ~exist('Titruth','var')
    theta=0:2*pi/NNodes:2*pi*(NNodes-1)/NNodes;
    Titruth=6*[cos(theta); sin(theta)];
end

t_node.Titruth=Titruth;
t_node=bearingNetworkAddMeasurements(t_node,'memberNode','Titruth','memberBearing','Yijtruth');
if isfield(t_node,'Er')
    t_node=bearingNetworkAddMeasurements(t_node,'memberEdges','Er','memberNode','Titruth','memberBearing','Yrijtruth');
end
