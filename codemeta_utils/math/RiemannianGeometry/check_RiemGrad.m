function check_RiemGrad(lf,f,gradf,x)
disp('Provided Gradient / Approximate Gradient')
disp([rot_vee(x,gradf(x)) rot_vee(x,approx_RiemGrad(lf,f,x))])
