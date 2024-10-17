Y0=stiefel_randn(eye(3,2))
Y0'*Y0
v0=stiefel_randTangentNormVector(Y0)
Y0'*v0
Yt=@(t) stiefel_exp(Y0,t*v0)
funPlot(@(t) Yt(t)'*Yt(t))