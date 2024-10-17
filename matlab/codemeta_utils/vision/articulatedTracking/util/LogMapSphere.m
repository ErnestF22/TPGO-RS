function mlog = LogMapSphere(x,y)

tol = 10^(-8);

if abs(norm(x,'fro')-1)>tol
     abs(norm(x,'fro')-1)
    error('Not a unit vector')
end

if abs(norm(y,'fro')-1)>tol
    abs(norm(y,'fro')-1)
    error('Not a unit vector')
end

x = x/norm(x);
y = y/norm(y);

trxy = x'*y;

geodist = acos(trxy) ;

if geodist< eps
    mlog = zeros(size(x));
else
    mlog = ((y-x*trxy)/sqrt(1-trxy^2)) *geodist;
  % mlog =  TODO WRITE THE LOGARITHM BETTER CONDITIONED
end

%mlog = (eye(size(x,1))-(x*x'))*mlog;

end

