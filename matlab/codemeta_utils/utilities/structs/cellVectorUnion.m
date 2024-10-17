function c=cellVectorUnion(a,b)
if ~iscell(a) || ~iscell(b)
    error('Inputs must be cell arrays')
end
la=length(a);
lb=length(b);
if la~=lb
    error('Inputs must be of the same length')
end
c=cell(size(a));
for ia=1:la
    c{ia}=union(a{ia},b{ia});
end
