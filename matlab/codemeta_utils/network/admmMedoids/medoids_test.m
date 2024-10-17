function medoids_test
resetRands()
nbPoints=100;
nbClusters=4;
dimPoints=2;
methodData='rand';
methodInit='zeros';
flagBias=false;
flagCenters=true;

optsMedoids={};

switch methodData
    case 'rand'
        x=rand(dimPoints,nbPoints);
    case 'uniform'
        x=repmat(linspace(0,1,nbPoints),dimPoints,1);
end
    
switch methodInit
    case 'rand'
        muInit=rand(dimPoints,nbClusters);
    case 'zeros'
        muInit=zeros(dimPoints,nbClusters);
end

if flagBias
    lambda=[repmat({[0;0]},1,nbClusters-2),[30;-30],[0;30]];
    optsMedoids=[optsMedoids {'bias',lambda}];
end

if flagCenters
    centers=0.1*[
        2 2 6 6;
        2 6 2 6
        ];
    centers=mat2cell(centers,2,ones(1,nbClusters));
    weights={{6},{0},{6},{20}};
    optsMedoids=[optsMedoids {'priorCenters',centers,weights}];
end    

[mu,output]=medoids(x,nbClusters,'muInit',muInit,'debug','verbose',optsMedoids{:});
itEnd=output.it;

figure(1)
plot(output.c)
switch dimPoints
    case 1
        figure(2)
        plot(itEnd,x,'bx');
        hold on
        plotPointsNumbers([itEnd*ones(size(mu));mu]);
        plot(0:itEnd,output.mu','r')
        hold off
    case 2
        figure(2)
        plotPoints(x,'bx');
        hold on
        plotPointsNumbers(mu);
        for iCluster=1:nbClusters
            plotPoints(squeeze(output.mu(:,iCluster,:)),'r-')
        end
        plotPoints(output.mu(:,:,end),'go')
        hold off
end

