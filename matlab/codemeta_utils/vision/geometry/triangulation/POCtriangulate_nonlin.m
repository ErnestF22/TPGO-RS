function POCtriangulate_nonlin
resetRands()

s=load('triangulate_test_dataset_data');
XTruth=s.X(:,1);
P=RTK2P(s.R,s.T,s.K);
xTruth=projectFromP(P,XTruth);

dX=randn(3,1);
X=@(t) XTruth+t*dX;

%check_der(@(t) errorVectDer(xTruth,P,X(t),dX), 'function')

X0=X(1);
opts=optimset('Display','Off','Jacobian','On');
[XEst,resnorm]=lsqnonlin(@(X) errorVect(xTruth,P,X), X0, [], [], opts);

disp(resnorm)
disp([X(0) XEst])

function [errVect,dErrVect]=errorVectDer(xTruth,P,X,dX)
[errVect,JErrVect]=errorVect(xTruth,P,X);
dErrVect=JErrVect*dX;

function [errVect,JErrVect]=errorVect(xTruth,P,X)
[err,JxProjected]=reprojectionError(xTruth,P,X);
errVect=err(:);
JErrVect=reshape(permute(squeeze(JxProjected),[1 3 2]),...
    size(JxProjected,1)*size(JxProjected,4),[]);
