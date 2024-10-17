%Normalize linear inequality constraints
%function [Cnorm,dnorm]=linearInequalitiesNormalize(C,d)
%Rescale the rows of C and d such that Cnorm*x>dnorm defines the same set
%as C*x>d, but each row of Cnorm has norm 1.
function [Cnorm,dnorm]=linearInequalitiesNormalize(C,d)
[Cnorm,l]=rnormalize(C);
dnorm=d./l;


