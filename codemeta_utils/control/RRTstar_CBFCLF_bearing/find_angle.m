function t = find_angle(a,b) %findLine
v = a-b;
if v(1)>=0 && v(2)>=0
    t = atand(v(2)/v(1));
end
if v(1)<0 && v(2)>0
    t = 180-atand(abs(v(2)/v(1)));
end
if v(1)<0 && v(2)<0
    t = 180+atand(v(2)/v(1));
end
if v(1)>0 && v(2)<0
    t = atand(v(2)/v(1));
end
end