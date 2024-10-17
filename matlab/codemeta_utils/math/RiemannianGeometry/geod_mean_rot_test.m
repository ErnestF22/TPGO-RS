R0=unif_random_rot(1);
v=randn(3,1);
allR=unif_random_rot(2);
[Ropt,topt]=geod_mean_rot(allR,R0,v)
