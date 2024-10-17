function [Res] = GF2dotproduct (V1,V2)
V = V1.*V2;
temp=rem(V,2);
Res=rem(sum(temp),2);
end