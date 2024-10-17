function POCrotateLocalization
dimData=3;
graphSampleNum=2;
nodeCentered=3;
[E,x,R]=locSegSampleNetworkGallery(graphSampleNum,dimData);
R0=rot_randn(eye(dimData));
NNodes=size(x,2);
x=x-x(:,nodeCentered)*ones(1,NNodes);
subplot(2,1,1)
graphShow(E,x,R)
axis equal
xNew=R0*x;
RNew=zeros(dimData,dimData,NNodes);
for iNode=1:NNodes
    RNew(:,:,iNode)=R0*R(:,:,iNode);
end
subplot(2,1,2)
graphShow(E,xNew,RNew)
axis equal
