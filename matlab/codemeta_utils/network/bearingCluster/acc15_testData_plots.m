function acc15_testData_plots

for nFigure=1:6
    switch nFigure
        case {1,2,3}
            load(['testData/testData' num2str(nFigure) '.mat'],'x','E')
        case {4,5,6}
            load('testData/testData4.mat','x')
            load(['testData/testData' num2str(nFigure) '.mat'],'E')
    end
    x(2,:)=max(x(2,:))-x(2,:);
    u=bearingCluster_getBearingsScalesFromE(x,E);
    membership=bearingCluster_clustering(E,u);

    %figure(nFigure)
    %bearingClustering_plot(x,E,membership)
    writeTikzFile(x,E,membership,nFigure)
    fprintf('testData%d: %d components\n',nFigure, length(unique(membership)))
end

function writeTikzFile(x,E,membership,nFigure)
fileDir='../figures/tikz';
fileNameBase='testData';

fileName=fullfile(fileDir,[fileNameBase num2str(nFigure) '.tex']);
fid=fopen(fileName,'wt');
bearingCluster_tikz(x,E,'fileId',fid,'flagNormalizeScale',true,'membership',membership);
fclose(fid);
