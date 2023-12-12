function uVec=bearingControlRotationVisibility(R,dx,y,nyMax,y0,funs)
dMax=1;
gradRPhiBound=bearingCostRotationVisibilityGradientRBoud(R,y,nyMax,y0,funs);
flagGradRSignAvailable=false;
if all(gradRPhiBound>=0)
    %gradRPhi is positive
    b=gradRPhiBound(1);
    flagGradRSignAvailable=true;
end
if all(gradRPhiBound<=0)
    %gradRPhi is negative
    b=gradRPhiBound(2);
    flagGradRSignAvailable=true;
end

if flagGradRSignAvailable
    if b==0
        dRVec=0;
    else
        dRVec=-dMax*b/b^2;
    end
    dx=[0;0];
else
    dRVec=0;
    g=bearingCostGeneral_gradient(y,R*y0,funs);
    dx=-dMax/(g'*g)*g;
end

uVec=[dRVec;dx];
