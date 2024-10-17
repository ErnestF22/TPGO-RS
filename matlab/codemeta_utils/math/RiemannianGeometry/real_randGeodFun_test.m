function real_randGeodFun_test
[xt,dxt,~,~,ddxt]=real_randGeodFun(randn(3,1),'quadraticSpeed');
check_der(xt,dxt)
check_der(dxt,ddxt)
