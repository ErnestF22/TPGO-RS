function T=G2T(G)
d=size(G,2)-1;
T=squeeze(G(1:d,d+1,:));
