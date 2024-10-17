function essential_logDiff_test
% resetRands()
switch 1
    case 1
        [Q1,~,~,~,v1Vec]=essential_randGeodFun();
        [Q2,~,~,~,v2Vec]=essential_randGeodFun();

        check_der(@(t) LogAndDer(Q1(t),Q2(t),v1Vec,v2Vec),'function',linspace(0,pi,100))
    case 2
        Q1=essential_randn();
        Q2=essential_randn();
        DLog=essential_logDiff(Q1,Q2);
        
        plotfuntrials(@(t) dLogRand(Q1,Q2,DLog),5)
end

function [Log,dLog]=LogAndDer(Q1,Q2,v1Vec,v2Vec)
[log,Q2r]=essential_log(Q1,Q2);
Log=essential_vee(Q1,log);
DLog=essential_logDiff(Q1,Q2r,'closestRepresentative');
dLog=DLog*[v1Vec;v2Vec];
 
% function [log,dlog]=logAndDer(Q1,Q2,v1Vec,v2Vec)
% [log,Q2r]=essential_log(Q1,Q2);
% Log=essential_vee(Q1,log);
% DLog=essential_logDiff(Q1,Q2r,'closestRepresentative');
% dLog=DLog*[v1Vec;v2Vec];
% dlog=essential_hat(Q1,dLog);

function v=dLogRand(Q1,Q2,DLog)
v1=essential_vee(Q1,essential_randTangentNormVector(Q1));
v2=essential_vee(Q2,essential_randTangentNormVector(Q2));
v=essential_tangentProjVerticalScalar(Q1,essential_hat(Q1,DLog*[v1;v2]));
