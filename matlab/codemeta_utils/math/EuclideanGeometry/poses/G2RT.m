function [R,T]=G2RT(G)
d=size(G,2)-1;
R=G(1:d,1:d,:);
T=squeeze(G(1:d,d+1,:));
