function POCConsistencyMAD
N=1e7;
x=sqrt(abs(randn(N,1)));
mu=median(x);
md=mad(x,1);
sd=std(x);
disp([mean(x) mu])
disp([sd md 1.374*md])
fprintf('%.8f\n',sd/md)
disp(sum(x<mu+2*md)/N)
