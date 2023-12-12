% check if the outputs of layers are met the constraint 
function y3=linearRestrict(x,A1,b1,A2,b2,A3,b3,s1,s2,s3)
[y1,y2,y3]=netForward(x,A1,b1,A2,b2,A3,b3);
if ~all(sign([y1';y2';y3])==[s1;s2;s3])
    y3=NaN(size(y3));
end
end