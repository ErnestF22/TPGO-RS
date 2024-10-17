function cumDistPerc(data,varargin)
data=shiftdim(data);
n=size(data,1);

data=sort(data);
plot(data,(1:n)/n*100,varargin{:})


