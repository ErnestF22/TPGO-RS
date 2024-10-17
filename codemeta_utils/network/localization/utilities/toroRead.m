%Read pose graph data from TORO into a testNetwork-compatible structure
%function t_node=toroRead(fileName)
function t_node=toroRead(fileName)
fid=fopen(fileName,'rt');
if fid<0
    error('Could not open file %s for reading')
end
data=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

NRecords=length(data{1});
NNodes=sum(strcmpi('VERTEX2',data{1}));
NEdges=sum(strcmpi('EDGE2',data{1}));

nodes.id=zeros(NNodes,1);
nodes.Ri=zeros(3,3,NNodes);
nodes.Ti=zeros(3,NNodes);

edges.E=zeros(NEdges,2);
edges.Rij=zeros(3,3,NEdges);
edges.Tij=zeros(3,NEdges);

nodesCnt=1;
edgesCnt=1;
for iRecord=1:NRecords
    switch data{1}{iRecord}
        case 'VERTEX2'
            nodes.id(nodesCnt)=data{2}(iRecord);
            nodes.Ti(1:2,nodesCnt)=[data{3}(iRecord);data{4}(iRecord)];
            nodes.Ri(:,:,nodesCnt)=theta2R(data{5}(iRecord));
            nodesCnt=nodesCnt+1;
        case 'EDGE2'
            edges.E(edgesCnt,:)=ids2idx(nodes.id,data{2}(iRecord),data{3}(iRecord));
            edges.Tij(1:2,edgesCnt)=[data{4}(iRecord);data{5}(iRecord)];
            edges.Rij(:,:,edgesCnt)=theta2R(data{6}(iRecord));
            edgesCnt=edgesCnt+1;
    end
end

t_node.NNodes=NNodes;
t_node.NEdges=NEdges;
t_node.E=edges.E;
t_node=testNetworkAddMeasurements(t_node,'method','given',RT2G(edges.Rij,edges.Tij));
t_node.gi=RT2G(nodes.Ri,nodes.Ti);

function R=theta2R(theta)
R=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

function e=ids2idx(allIds,id1,id2)
e=[find(allIds==id1),find(allIds==id2)];
if length(e)~=2
    error('Error on an edge of the TORO file')
end
