function essential_fromE_test
for seed=0:10
    fprintf('Seed: %d\n',seed)
    resetRands()
    Q=essential_randn();
    E=essential_getE(Q);
    Q1=essential_fromE(E);
    E1=essential_getE(Q1);

    disp([Q Q1])
    disp(essential_dist(Q,Q1))
    disp([E E1])
end
