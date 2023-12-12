% forward network function for NN with 3 layers, the output is scalar
function y = netForward(Ab_set,z,x)
L = size(Ab_set,2);
[y1, z1] = forward(x,L,Ab_set);
if ~all(z==z1)
    y = NaN(size(y1));
else
    y = y1;
end
end