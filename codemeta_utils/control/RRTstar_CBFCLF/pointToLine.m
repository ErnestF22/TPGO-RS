function dis = pointToLine(obs,n,p)
t1 = p-obs;
t2 = p-n;
t = (t1'*t2)/(norm(t2,2)^2);
if ((t<=1) && (t>=0))
    d0 = sum(abs(cross([t1(1) t1(2) 0],[t2(1) t2(2) 0])));
    dis = d0/norm(t2,2);
else
    dis = Inf;
end
end