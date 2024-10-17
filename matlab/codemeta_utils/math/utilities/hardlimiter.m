%function x=hardlimiter(x,t1,t2)
%Clips the value of x to be between t1 and t2. If omitted, t1=1 and t2=-t1

%%AUTORIGHTS%%

function x=hardlimiter(x,t1,t2)
if exist('t1','var')==0
    t1=1;
end

if exist('t2','var')==0
    t2=-t1;
end

if x>t1
    x=t1;
elseif x<t2;
    x=t2;
end
