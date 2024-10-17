function d=hyperbolic_lorentzMetric(v1,v2)
d=v1(1:end-1)'*v2(1:end-1)-v1(end)*v2(end);


