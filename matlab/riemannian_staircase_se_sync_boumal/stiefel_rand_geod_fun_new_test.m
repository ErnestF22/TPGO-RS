function stiefel_rand_geod_fun_new_test()
Y0=stiefel_randn(eye(3,2));
[Yt,Ht]=stiefel_rand_geod_fun_new(Y0);
figure(1)
funPlot(@(t) Yt(t)'*Yt(t)-eye(2))
title('Should have all lines near zero')
figure(2)
funPlot(@(t) diag(Ht(t)'*Yt(t)))
title('Should have all lines near zero')
figure(3)
funPlot(@(t) trace([0 1;1 0]'*Ht(t)'*Yt(t)))
title('Should have all lines near zero')
funCheckDer(Yt,Ht)
title('Numerical check of the tangent')