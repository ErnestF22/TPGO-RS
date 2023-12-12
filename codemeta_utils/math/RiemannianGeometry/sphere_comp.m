%Returns the equivalent of x1 o x2, the composition of two elements of the
%group. If we pass the optional argument 'inv', compute x1^-1 o x2
function [x12]=sphere_comp(x1,x2,varargin)
flaginv=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case {'inv','inverse'}
            flaginv=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

D=size(x1,1);
e=[1;zeros(D-1,1)];
h1=sphere_log(e,x1);
if(flaginv)
    h1=-h1;
end
h2=sphere_log(e,x2);
h1at2=sphere_parallel(e,h2,h1);
x12=sphere_exp(x2,h1at2);
