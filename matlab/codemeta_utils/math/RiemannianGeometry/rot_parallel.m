%function [d1]=rot_parallel(R,h,d,varargin)
%Parallel transport d, a vector in the tangent space at y, along the
%geodesic given by h
%If the optional argument 'toRotation' is passed, h is considered to be
%another rotation and the geodesic along which the function parallel
%transports is identified by the endpoints Rstart=R and Rend=h
function [d1]=rot_parallel(R,h,d,varargin)
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'torotation'
            h=rot_log(R,h);             %h is not a vector but the destination rotation. Compute the opportune h
    end
    ivarargin=ivarargin+1;
end

%for the method used see Edelman et al., eq (2.18)
b=R'*d;     %obtain skew sym representation
h=R'*h;     %ditto for h
d1=R*rot_exp(eye(3),h/2)*b*rot_exp(eye(3),h/2);


