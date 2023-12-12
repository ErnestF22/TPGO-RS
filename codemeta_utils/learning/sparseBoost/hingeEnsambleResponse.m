%Evaluate reponse of a linear ensamble of hinges
%function r=hingeEnsambleResponse(x,w,v,a,b)
function r=hingeEnsambleResponse(x,w,v,a,b)
[d,Nw]=size(w);
flagAppendOpposite=length(a)==2*Nw;
r=hingeResponse(x,w,v,'flagAppendOpposite',flagAppendOpposite);
r=r*a+b;
