%Plots an image and overlaps points on it
%function [hImage,hPlot]=sfm_rawDisplayFeature(img,x,varargin)
%Inputs
%   img     array with the image data or string with the file name
function [hImage,hPlot]=sfm_rawDisplayFeature(img,x,varargin)
optsPlot={};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optsplot'
            ivarargin=ivarargin+1;
            optsPlot=[optsPlot{:} varargin{ivarargin}];
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ischar(img)
    img=imread(img);
end

h(1)=imshow(img);
flagIsHold=ishold();
hold on
h(2)=plot(x(1,:),x(2,:),'*','Color',[0 0.8 0], optsPlot{:});
if ~flagIsHold
    hold off
end

if nargout>0
    hImage=h(1);
    if nargout>1
        hPlot=h(2);
    end
end
