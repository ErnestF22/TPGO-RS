function p=intersectionTwoLines(a1,b1,a2,b2)
p = [0;0];
p(1) = (b2-b1)/(a1-a2);
p(2) = p(1)*a1+b1;
if (p(1)==Inf) || (p(1)==-Inf)
    p = [1000;1000];
end
end