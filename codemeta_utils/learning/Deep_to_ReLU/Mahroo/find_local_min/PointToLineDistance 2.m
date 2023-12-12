function d = PointToLineDistance(a,b,c,p)
% ax+by+c=0
d = abs(a*p(1)+b*p(2)+c)/(sqrt(a^2+b^2+1));
end