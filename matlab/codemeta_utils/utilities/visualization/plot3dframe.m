%Plots 3-D axes
%function plot3dframe(o,r,style)
%Inputs
%   o   column vector for the origin of the frame
%   r   the column of R are the axis vectors (they should be orthonormal)

%%AUTORIGHTS%%

function plot3dframe(o,r,style,varargin)
labels={'x','y','z'};
scale=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
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

r=max(scale)*r;

h=ishold();
if(exist('style','var') && ~isempty(style))
    plot3dvect(o,r(:,1),labels{1},style)
    hold on
    plot3dvect(o,r(:,2),labels{2},style)
    plot3dvect(o,r(:,3),labels{3},style)
else
    plot3dvect(o,r(:,1),labels{1},'r')
    hold on
    plot3dvect(o,r(:,2),labels{2},'g')
    plot3dvect(o,r(:,3),labels{3},'b')
end
if(~h)
    hold off
end

