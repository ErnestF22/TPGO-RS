function kmeans_plot(xi,mu,idx,muInit,lambda,z)
if ~exist('idx','var') || isempty(idx)
    idx=kmeans_augmented2_idxGroup(xi,mu);
end
plotGroups(xi,idx)
hold on
plotPoints(mu,'x')
plotPoints(muInit,'kx')
for iMean=1:size(mu,2)
    text(mu(1,iMean),mu(2,iMean),['\mu_' num2str(iMean)])
    if ~isempty(muInit)
        text(muInit(1,iMean),muInit(2,iMean),['\mu_{0' num2str(iMean) '}'])
    end
end
if exist('z','var') && ~isempty(z)
    plotPoints(z,'d')
    for direction=1:2
        for iMean=1:size(z,3)
            text(z(1,direction,iMean),z(2,direction,iMean),['z_' num2str(iMean)])
        end
        quiver(mu(1,:),mu(2,:),squeeze(lambda(1,direction,:))',squeeze(lambda(2,direction,:))')
    end
    hold off
end
