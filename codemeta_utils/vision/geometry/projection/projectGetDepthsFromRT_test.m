function projectGetDepthsFromRT_test
load('triangulate_test_dataset_datacalibrated')
XTransformed=rigidTransform(R,T,X,'poses');
lambda=projectGetDepthsFromRT(R,T,X,'poses');
x=projectFromRT(R,T,X,'poses');
xHom=homogeneous(x,3);
XTransformedBackprojected=repmat(permute(lambda,[3 2 1]),3,1).*xHom;

disp([XTransformed;XTransformedBackprojected])
