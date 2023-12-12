function realDyn_statePack_test
NT=5;
T=randn(3,NT);
v=randn(3,NT);
[T1,v1]=realDyn_stateUnpack(realDyn_statePack(T,v));
disp(norm(T(:)-T1(:)))
disp(norm(v(:)-v1(:)))

