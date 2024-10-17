function v=sphere_randTangentNormVector(y,varargin)
flagPermute=false;
if size(y,3)==1
    y=permute(y,[1 3 2]);
    flagPermute=true;
end
v=lie_randTangentNormVector(sphere_funs(),y,varargin{:});
if flagPermute
    v=permute(v,[1 3 2]);
end
