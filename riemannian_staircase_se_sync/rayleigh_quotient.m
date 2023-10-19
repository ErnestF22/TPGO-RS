function retval = rayleigh_quotient(x, A)
%RAYLEIGH_QUOTIENT Returns Rayleigh Quotient computed between vector x and
% matrix A 
% i.e., 
% retval = (x' * A * x) / (x' * x)

%TODO: insert correct dimension checks (and possibly others)

retval = (x' * A * x) / (x' * x);

end

