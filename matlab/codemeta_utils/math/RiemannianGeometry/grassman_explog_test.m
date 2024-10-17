function grassman_explog_test
n=7;
p=3;

Y1=grassman_randn(eye(7,3),1);
Y2=grassman_randn(eye(7,3),1);

disp('iterative')
disp(subspace(Y2,grassman_exp(Y1,grassman_log(Y1,Y2,'method','iterative'))))

disp('closed form')
disp(subspace(Y2,grassman_exp(Y1,grassman_log(Y1,Y2,'method','gsvd'))))
