%Refine estimation of essential matrix by 
function E=epipolarEstimateESampson(x1,x2,EInit)

QInit=essential_fromE(EInit);
Q=essential_minimizeEpipolarCostSampsonSq(cat(3,x1,x2),QInit);
E=essential_toE(Q);
