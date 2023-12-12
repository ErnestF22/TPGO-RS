%Pass from bearing vectors to cosines
%function c=bearingComputeCosine(y,yg)
%Inputs
%   y    current measured bearing vectors
%   yg   bearing vectors at the goal position
function c=bearingComputeCosine(y,yg)
NY=size(y,2);
NYg=size(yg,2);
if NYg==1 && NYg<NY
    yg=yg*ones(1,NY);
end
if NY==1 && NY<NYg
    y=y*ones(1,NYg);
end

c=min(1,max(-1,sum(yg.*y)));
