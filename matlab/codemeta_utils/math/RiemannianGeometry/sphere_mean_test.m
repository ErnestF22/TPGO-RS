N=1000;
D=7;
sigmanoise=0.1;

e=[1;zeros(D-1,1)];

y=e*ones(1,N)+sigmanoise*randn(D,N);

y=y./(ones(D,1)*sqrt(sum(y.^2)));

ymean=sphere_mean(y);

display(ymean)

