function hyperbolic_test
resetRands();

n=3;
N=5;

x=[zeros(n,1); 1];

v=hyperbolic_tangentProj(x,randn(n+1,N));

y=hyperbolic_exp(x,v);

display(v)
display(y)

disp('v-v reconstr')
disp([v-hyperbolic_log(x,hyperbolic_exp(x,v))])
disp('y-y reconstr')
disp([y-hyperbolic_exp(x,hyperbolic_log(x,y))])

disp('distances of y from y(:,:,1)')
disp(hyperbolic_dist(y(:,1),y))

ymean=hyperbolic_mean(permute(y,[1 3 2]),'showcost','showit');
display(ymean)