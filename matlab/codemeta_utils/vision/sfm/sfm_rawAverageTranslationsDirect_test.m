function sfm_rawAverageTranslationsDirect_test
t_node=testNetworkBuildTestNetwork();
E=testNetworkGetEdges(t_node);
Ri=t_node.Ritruth;
Tij=t_node.Tijtruth;

TiRef=t_node.Titruth;
TiRef=TiRef-mean(TiRef,2)*ones(1,size(TiRef,2));
TiRef=TiRef/norm(TiRef(:));
save temp

Ti=sfm_rawAverageTranslationsDirect(Ri,Tij,E);
disp([Ti;TiRef])
disp(norm(Ti-TiRef,'fro'))
