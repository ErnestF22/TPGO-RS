function l=translationLogLikelihood(Ri,Ti,Tj,Tij,Gamma)
v=Tij-Ri'*(Tj-Ti);
l=(v'*Gamma*v)/2;
