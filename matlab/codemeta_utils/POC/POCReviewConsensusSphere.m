function POCReviewConsensusSphere
x=sphere_randn();
s=randn(3,1);

disp([norm(cross(x,s))^2 s'*s-(s'*x)^2]);
