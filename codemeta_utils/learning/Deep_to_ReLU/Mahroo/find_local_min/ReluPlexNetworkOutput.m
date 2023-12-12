%Given the input, and a set of weigthed matrices and bias matrices, this
%function returns which region input belongs to. 

function [y,layer] = ReluPlexNetworkOutput(x,layer,L)
% x is the input to the network and y is the output of the network
[y,layer]=ReluPlexNetworkForward(x,layer,L);
end