function admm_averagingKmeans_plot(nodes,xi0,allMu)
d=size(xi0{1},1);
xi0Group=cell2groupIdx(xi0);
plotGroups([xi0{:}],[xi0Group{:}])
hold on
nbNodes=length(nodes);
K=size(nodes(1).mu,3);
for iNode=1:nbNodes
    for k=1:K
        mu=nodes(iNode).mu(:,1,k);
        plotPoints(mu,'x')
        text(mu(1),mu(2),['\mu:i' num2str(iNode) 'k' num2str(k)])
        zINode=nodes(iNode).z;
        for j=1:size(zINode,2)
            z=zINode(:,j,k);
            plotPoints(z,'d')
            text(z(1),z(2),['z:i' num2str(iNode) 'j' num2str(j) 'k' num2str(k)])
            lambdaINode=nodes(iNode).lambda;
            for direction=1:2
                l=lambdaINode(:,j,k);
                quiver(z(1),z(2),l(1),l(2),0.1)
            end
        end
    end
end
if exist('allMu','var')
    switch d
        case 2
            plot(squeeze(allMu(:,1,:)),squeeze(allMu(:,2,:)))
        case 3
            plot3(squeeze(allMu(:,1,:)),squeeze(allMu(:,2,:)),squeeze(allMu(:,3,:)))
        otherwise
            error('Cannot display this number of dimensions')
    end
    mu0=squeeze(allMu(1,:,:,:));
    plotPoints(mu0,'kx')
end

hold off
axis equal
