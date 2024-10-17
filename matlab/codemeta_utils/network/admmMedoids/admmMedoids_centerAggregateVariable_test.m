function admmMedoids_centerAggregateVariable_test
problem=admmMedoids_dataset();
nbNodes=problem.nbNodes;
nbClusters=problem.nbClusters;
A=adjgallery(nbNodes,'banded',2);
nodeData=admmMedoids_initNodeData(problem,A);
nbDim=size(nodeData(1).x,1);
for iNode=1:nbNodes
    for idxJNode=1:length(nodeData(iNode).idxNeighbors)
        jNode=nodeData(iNode).idxNeighbors(idxJNode);
        for k=1:nbClusters
            for d=1:nbDim
                for l=1:2
                    nodeData(iNode).lambda(d,idxJNode,k,l)=sum([iNode jNode k l d].*(10.^(4:-1:0)));
                end
                nodeData(iNode).z(d,idxJNode,k)=sum([iNode jNode k d].*(10.^(3:-1:0)));
            end
        end
    end
end

for iNode=1:2
    disp(['Node ' num2str(iNode) ' neighbors'])
    disp(nodeData(iNode).idxNeighbors)
    disp(['Node ' num2str(iNode) ' lambda'])
    disp(admmMedoids_centerAggregateVariable(nodeData,iNode,'lambda'))
    disp(['Node ' num2str(iNode) ' zeta'])
    disp(admmMedoids_centerAggregateVariable(nodeData,iNode,'z'))
end

fprintf(['Meaning of each digit in lambda\n',...
    '\tStart node\n'...
    '\tEnd node\n'...
    '\tCluster\n'...
    '\tEndpoint\n'...
    '\tDimension\n'])
disp('For zeta it is the same without the endpoint digit')

disp('Checks:')
fprintf('Either first or second digit of lambda and zeta is the same as the node number: ')
for iNode=1:nbNodes
    data=vec(admmMedoids_centerAggregateVariable(nodeData,iNode,'lambda'))';
    for d=data
        assert(get_digit(d,4)==iNode || get_digit(d,5)==iNode)
    end
    data=vec(admmMedoids_centerAggregateVariable(nodeData,iNode,'z'))';
    for d=data
        assert(get_digit(d,3)==iNode || get_digit(d,4)==iNode)
    end
end
disp('Ok')

fprintf('Second dimension of gathered lambda and zeta is twice the number of neighbors: ')
for iNode=1:nbNodes
    data=admmMedoids_centerAggregateVariable(nodeData,iNode,'lambda');
    assert(size(data,2)==2*length(nodeData(iNode).idxNeighbors))
end
disp('Ok')

fprintf('Third digit of lambda and zeta is the same as the cluster number: ')
for iNode=1:nbNodes
    data=admmMedoids_centerAggregateVariable(nodeData,iNode,'lambda');
    for iCluster=1:size(data,3)
        assert(all(vec(get_digit(data(:,:,iCluster,:),3)==iCluster)))
    end
    data=admmMedoids_centerAggregateVariable(nodeData,iNode,'z');
    for iCluster=1:size(data,3)
        assert(all(vec(get_digit(data(:,:,iCluster),2)==iCluster)))
    end
end
disp('Ok')


% fprintf(['Verify that:\na) values for z are always 10^4 more than values for lambda\n' ...
%     'b) either value on the fourth or third digit from the right is equal to the iNode number\n' ...
%     'c) third digit from the left is equal to the cluster number\n'])

%Function to return the exponent-th digit from the left of the decimal
%point
function x_digit=get_digit(x,exponent)
mask=10^(exponent-1);
x_digit=mod(floor(x/mask),10);
