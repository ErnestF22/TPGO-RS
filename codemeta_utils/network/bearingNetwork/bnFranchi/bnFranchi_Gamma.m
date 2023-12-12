function gammalmn=bnFranchi_Gamma(bmn,bml,blm,bln)
sln=sqrt(1-bearingComputeCosine(bmn,bml).^2);
sml=sqrt(1-bearingComputeCosine(bln,bnFranchi_RzFromBearings(blm,bml)*bmn).^2);
gammalmn=sln/sml;

