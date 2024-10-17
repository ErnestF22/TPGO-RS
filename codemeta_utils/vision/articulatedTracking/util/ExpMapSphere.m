function mexp = ExpMapSphere(x,xi)

tol = 10^(-5);

if abs(norm(x)-1)>tol
  %  abs(norm(x)-1)
  %  error('Not a unit vector')
end

nxi = norm(xi(:));

if nxi < eps
    mexp = x;
else
    mexp = x*cos(nxi) + xi*(sin(nxi)/nxi);
end

mexp = mexp / norm(mexp);

end

