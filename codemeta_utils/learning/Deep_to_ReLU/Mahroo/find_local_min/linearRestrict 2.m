% check if the outputs of layers are met the constraint 
function y=linearRestrict(x,layer,L,S)
[y,layer]=ReluPlexNetworkForward(x,layer,L);
Y=[];
for i=1:L
    Y = [Y;layer(i).y];
end
if ~all(sign(Y)==S)
    y=NaN(size(y));
end
end