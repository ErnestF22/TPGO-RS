function ab=multiprodN(a,b,varargin)
ab=multiprod(a,b);
for ivarargin=1:length(varargin)
    ab=multiprod(ab,varargin{ivarargin});
end
end %file function