% compute the output of network and corresponding z
function [y, z] = forward(x,L,Ab_set)
z = [];
for i =1:L
    h = hidden_layer(Ab_set(i).A,Ab_set(i).b,x);
    
    hr = round(h,1);
    z = [z; (hr>=0)];
    
    yi = Re(h);
    
    x = yi;
end
y = yi;
end