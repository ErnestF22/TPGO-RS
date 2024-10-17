%Plots 2-D axes
%function plot2dframe(o,R,style,varargin)
%Inputs
%   o  column vector for the origin of the frame
%   R  the column of R are the axis vectors (they should be orthonormal)
function plot2dframe(o,r,style,varargin)
labels={'x','y','z'};
scale=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'scale'
            ivarargin=ivarargin+1;
            scale=varargin{ivarargin};
        case 'nolabels'
            labels={'','',''};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

r=scale*r;

h=ishold();
if(exist('style','var') && ~isempty(style))
    plot2dvect(o,r(:,1),labels{1},style)
    hold on
    plot2dvect(o,r(:,2),labels{2},style)
else
    plot2dvect(o,r(:,1),labels{1},'r')
    hold on
    plot2dvect(o,r(:,2),labels{2},'g')
end
if(~h)
    hold off
end
