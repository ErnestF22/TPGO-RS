function triang_exp_test
X=[1;1;2];
x=triang_projectX(X);

switch 2
    case 1
        v=randn(2,2);
        v=triang_tangProj(x,v);

        disp([v triang_log(x,triang_exp(x,v))])
    case 2
        T=triang_tangentBasis(x);
        f=@() triang_epipolarConstraint(triang_exp(x,reshape(T*randn(3,1),2,2)));
        plotfuntrials(f)
end

