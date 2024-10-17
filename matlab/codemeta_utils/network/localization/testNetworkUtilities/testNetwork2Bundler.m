function br=testNetwork2Bundler(t_node,X,dirName)
structType=testNetworkDetectStructType(t_node);

NCameras=testNetworkGetNumberOfNodes(t_node);
NPoints=size(X,2);
br.dirName=dirName;
br.fileHeader='# Bundle file v.03';
br.NCameras=NCameras;
br.NPoints=NPoints;

cm.f=ones(NCameras,1);
cm.k1=zeros(NCameras,1);
cm.k2=zeros(NCameras,2);

[cm.R, cm.T]=testNetworkGetRotTransl(t_node,'fieldname','gitruth');
cm.X=cell(1,NCameras);
[cm.X{:}]=deal(X);

switch structType
    case 'array'
        for iNode=1:length(t_node)
            cm.xy=cell(1,NCameras);
            cm.xy{iNode}=t_node(iNode).ximage;
        end
    case 'single'
        cm.xy=t_node.ximage;
end

br.cameras=cm;
pt.position=X;
pt.color=zeros(size(X));

vl=cell(1,NPoints);
for iPoint=1:NPoints
    view.camera=1:NCameras;
    view.idx=iPoint*ones(1,NCameras);
    view.xy=zeros(2,NCameras);
    for iCamera=1:NCameras
        view.xy(:,iCamera)=cm.xy{iCamera}(1:2,iPoint);
    end
    vl{iPoint}=view;
end
pt.NViews=NCameras*ones(1,NPoints);

pt.viewList=vl;
br.points=pt;
