function [y,layer]=ReluPlexNetworkForward(x,layer,L)
for i =1:L
    h = hidden_layer(layer(i).A,layer(i).b,x);
    layer(i).z = (h>0);
    layer(i).y = ReLU(h);
    x = layer(i).y;
end
y = x;
end