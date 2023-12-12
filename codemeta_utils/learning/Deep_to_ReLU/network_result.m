function y3 = network_result(A1,b1,A2,b2,A3,b3,s1,s2,s3,x)
[z1,z2,z3,y3]=forward_net(A1,b1,A2,b2,A3,b3,x);
if ~all([z1;z2;z3]==[s1;s2;s3])
    y3=NaN(size(y3));
else
    y3=y3;
end
end