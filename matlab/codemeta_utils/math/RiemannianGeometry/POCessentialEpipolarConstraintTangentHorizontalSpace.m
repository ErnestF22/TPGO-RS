function POCessentialEpipolarConstraintTangentHorizontalSpace
[G,X]=testNetworkCreateAbsolutePoses(3);
G=G(:,:,1:2);
x=projectFromG(G,X);
Q0=essential_fromG(G(:,:,1),G(:,:,2),'poses');
Qt=essential_randGeodFun(Q0);

plotfun(@(t) gradVertComp(Qt,x,t),'angle')

function a=gradVertComp(Qt,x,t)
Q=Qt(t);
[~,De]=essential_evaluateEpipolarCost(Q,x(:,:,1),x(:,:,2));
g=sum(De,2);
a=essential_tangentProjVerticalScalar(Q,essential_hat(Q,g));
