function [z1,z2,z3,y3]=forward_net(A1,b1,A2,b2,A3,b3,x)
h1 = hidden_layer(A1,b1,x);
y1 = ReLU(h1);
z1 = int8((y1>0)');
h2 = hidden_layer(A2,b2,y1');
y2 = ReLU(h2);
z2 = int8((y2>0)');
h3 = hidden_layer(A3,b3,y2');
y3 = ReLU(h3);
z3 = int8((y3>0));
end