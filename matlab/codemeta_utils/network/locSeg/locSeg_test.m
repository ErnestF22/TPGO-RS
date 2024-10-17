function locSeg_test(graphSampleNum,constraintSampleNum,dimData)
resetRands()
if ~exist('dimData','var')
    dimData=3;
end
if ~exist('graphSampleNum','var')
    graphSampleNum=9;
end
if ~exist('constraintSampleNum','var')
    constraintSampleNum=3;
end

[E,x,R,constraintList,constraintParametersList]=locSegSampleNetworkGallery(graphSampleNum,'dimData',dimData,'constraintNum',constraintSampleNum);
figure(1)
graphShow(E,x,R)
%[constraintList,constraintParametersList,E]=locSegSampleConstraintsGallery(constraintSampleNum,E,x,R);
locSeg(x,R,E,constraintList,constraintParametersList)
