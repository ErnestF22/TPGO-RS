function [Res] = GF2sum (V1,V2)
V = V1+V2;
Res=rem(V,2);
end