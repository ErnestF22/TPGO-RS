%Run a function on clusters of data
%function r=indexedFun(m,x,fun)
%Split the vector x according to the membership indeces in m and run f on
%each cluster
function [mUnique,r]=indexedFun(m,x,fun)
mUnique=unique(m);
r=zeros(size(mUnique));
for im=1:length(mUnique)
    r(im)=fun(x(m==mUnique(im)));
end