%function g=geom_mean(x)
%Compute the geometric mean of a vector x
function g=geom_mean(x,w)
if(nargin<2)
    w=ones(size(x));
end
g=(prod(x.^w))^(1/sum(w));
