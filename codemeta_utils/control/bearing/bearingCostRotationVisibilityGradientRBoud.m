function gradRPhiVecBound=bearingCostRotationVisibilityGradientRBoud(R,y,nyMax,y0,funs)
a=bearingCostRotationVisibilityGradientRTerms(R,y,y0,funs);
a=a*nyMax;
gradRPhiVecBound(1)=zeroIfEmpty(sum(a(a>0)));
gradRPhiVecBound(2)=zeroIfEmpty(sum(a(a<0)));

%Return zero if x is empty, otherwise return x itself
function x=zeroIfEmpty(x)
if isempty(x)
    x=0;
end
