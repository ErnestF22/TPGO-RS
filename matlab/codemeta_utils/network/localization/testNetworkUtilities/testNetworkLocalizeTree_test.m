function testNetworkLocalizeTree_test
randn('state',0)
rand('state',0)

flagDisplayResult=false;

testNum=1;
switch testNum
    case 1
        methodAbsolutePoses='references';
        structType='array';
    case 2
        methodAbsolutePoses='poses';
        structType='array';
    case 3
        methodAbsolutePoses='references';
        structType='single';
    case 4
        methodAbsolutePoses='poses';
        structType='single';
    otherwise
        error();
end

N=1000;
A=adjgallery(N,'banded',30);
[G,X]=testNetworkCreateAbsolutePoses(N);

t_node=testNetworkCreateStruct(A,'Type',structType);  
t_node=testNetworkAddGroundTruth(t_node,G,methodAbsolutePoses);
t_node=testNetworkAddMeasurements(t_node,'Method','Truth');
%t_node=testNetworkAddMeasurements(t_node,'Method','Truth','Unnormalized');
%t_node=testNetworkInitializeStates(t_node,'MethodR','Truth','MethodT','Truth','MethodScale','Truth');

ETree=testNetworkGenerateRandomTree(t_node);
t_node=testNetworkLocalizeTree(t_node,ETree,methodAbsolutePoses);

t_node=testNetworkCompensate(t_node,methodAbsolutePoses);

if flagDisplayResult
    figure(1)
    testNetworkDisplay(t_node,'DisplayEdges',methodAbsolutePoses);
    hold on
    testNetworkDisplay(t_node,'Member','gi','Estimated',methodAbsolutePoses);
    hold off
end

[rotErr,translErr,scaleRatios]=testNetworkComputeErrors(t_node,methodAbsolutePoses);
[rotErrMs,translErrMs,scaleRatiosMs]=testNetworkComputeErrors(t_node,'Measurements',methodAbsolutePoses);
