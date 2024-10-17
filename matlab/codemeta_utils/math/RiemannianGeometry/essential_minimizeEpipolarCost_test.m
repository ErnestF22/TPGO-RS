function essential_minimizeEpipolarCost_test
%resetRands(1)
NRestarts=5;
NX=5;
[G,X]=testNetworkCreateAbsolutePoses(3);
X=X(:,1:NX);
G=G(:,:,1:2);
%testNetworkDisplay(G,'points',X);
x=projectFromG(G,X);
%plot(x(1,:,1),x(2,:,1),'*')
QTruth=essential_fromG(G(:,:,1),G(:,:,2),'poses');


fprintf('#')
c=NaN(1,NRestarts);
for nRestart=1:NRestarts
    fprintf('.')
    Q0=essential_randn();
    QEst=essential_minimizeEpipolarCost(Q0,x(:,:,1),x(:,:,2));
    c(nRestart)=sum(essential_evaluateEpipolarConstraint(QEst,x(:,:,1),x(:,:,2)).^2);
    if c(nRestart)<1e-9
        disp('Success')
        break;
    end
end
fprintf('\n')

disp('Min cost for trials')
disp(min(c(~isnan(c))))
disp('Distance between QTruth/Q0 and QTruth/QEst')
disp([essential_dist(QTruth,Q0) essential_dist(QTruth,QEst)])
disp('Epipolar constraint evaluation')
disp(essential_evaluateEpipolarConstraint(QEst,x(:,:,1),x(:,:,2)));
xHom=homogeneous(x);
ETruth=essential_getE(QTruth);
EEst=essential_getE(QEst);

save([mfilename '_data'])
