function [A,b] = fitLine(t,p)
if mod(t,90) == 0 %vertical line
    A = [1 0];
    b = p(1);
end
if mod(t,180) == 0 %horizontal line
    A = [0 1];
    b = p(2);
end

if mod(t,90) ~= 0
    slop = tand(t);
    A = [slop 1];
    b = p(2)-(slop)*p(1);
end
end


