%Same as ind2sub but produces a vector instead of separate variables
function out = myind2sub(siz,ndx)
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1,
  v = floor(ndx/k(i))+1;  
  out(i) = v; 
  ndx = rem(ndx,k(i));
end
end