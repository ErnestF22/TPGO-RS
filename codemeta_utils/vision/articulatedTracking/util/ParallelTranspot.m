function zetat = ParallelTranspot(x,xt,zeta0)

% Parallel transport sphere

if any(isnan(x)) | any(isnan(xt)) | any(isnan(zeta0))
    error('Something is not a number')
end

tol = 10^(-8);

if abs(norm(x,'fro')-1)>tol
    norm(x,'fro')
    error('Not a unit vector')
end

if abs(norm(xt,'fro')-1)>tol
    abs(norm(xt,'fro')-1)
    error('Not a unit vector')
end


zeta0 = (eye(3)-x*x')*zeta0;


% transports zeta0 from the tangent space of x to the tangent space of xt

zetat = zeta0;

xi =  LogMapSphere(x,xt);
xinorm = norm(xi);

if xinorm > eps
    u = xi/xinorm;
    zetat = -x*sin(xinorm)*(u'*zeta0) + u*cos(xinorm)*u'*zeta0 + (eye(size(x,1))-u*u')*zeta0;
end



end