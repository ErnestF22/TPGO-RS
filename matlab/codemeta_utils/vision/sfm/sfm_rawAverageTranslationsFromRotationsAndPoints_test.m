function sfm_rawAverageTranslationsFromRotationsAndPoints_test
nDataset=1;

switch nDataset
    case 1
        t_node=testNetworkBuildTestNetwork();
        E=testNetworkGetEdges(t_node);
        x1=t_node.ximage(E(:,1));
        x2=t_node.ximage(E(:,2));
        Ri=t_node.Ritruth;
        TiRef=t_node.Titruth;
    case 2
        load('sfm_test_data_2','data')

        E=(sfm_getMatchIdxImg(data,'matchFiltered'))';
        [x1,x2]=sfm_getFeatureLocationFromMatchId(data,[],'normalized');
        [Ri,TiRef]=G2RT(data.poseTruth);
end

disp('# Median/max algebraic error for each edge')
for iEdge=1:size(E,1);
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    EEdge=epipolarBuildEFromRT(Ri(:,:,iNode),TiRef(:,iNode),Ri(:,:,jNode),TiRef(:,jNode));
    err=epipolarResiduals(EEdge,x1{iEdge},x2{iEdge},'algebraicabs');
    fprintf('    %.4e  %.4e\n',median(err),max(err))
end

TiRef=TiRef-mean(TiRef,2)*ones(1,size(TiRef,2));
TiRef=TiRef/norm(TiRef(:));

[Ti,A]=sfm_rawAverageTranslationsFromRotationsAndPoints(Ri,x1,x2,E);
disp('# Results')
disp('[Ti;TiRef]')
disp([Ti;TiRef])
disp('Frobenious norm of error')
disp(norm(Ti-TiRef,'fro'))
cumDistPerc(abs(A*TiRef(:)))
hold on
cumDistPerc(abs(A*Ti(:)),'r')
hold off
