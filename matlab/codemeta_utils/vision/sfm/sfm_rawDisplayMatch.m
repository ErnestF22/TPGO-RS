%Show two images overlapped and corresponding points
%function sfm_rawDisplayMatch(img1,img2,x1,x2,varargin)
%By default, adds a slider to change the relative opacity of the two
%overlapping images.
%Inputs
%   img1,img2   images to display
%   x1,x2       two [2 x NPoints] arrays with coordinates of the points to
%               display
%Optional arguments
%   'alphaWeights',aw   [2 x 1] vector containing the relative transparencies
%                       of the two images (it is automatically normalized)
%   'noslider'          do not show the slider
function sfm_rawDisplayMatch(img1,img2,x1,x2,varargin)
alphaWeights=[0.5 0.5];
flagSlider=true;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'alphaweights'
            ivarargin=ivarargin+1;
            alphaWeights=varargin{ivarargin};
        case 'flagslider'
            ivarargin=ivarargin+1;
            flagSlider=varargin{ivarargin};
        case 'noslider'
            flagSlider=false;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

alphaWeights=alphaWeights/sum(alphaWeights);
cla
[hImg,hPlot]=sfm_rawDisplayFeature(img1,x1,'optsPlot',{'Color',[0.9 0 0]});
set(hImg,'AlphaData',alphaWeights(1));
hold on
hImg=sfm_rawDisplayFeature(img2,x2,'optsPlot',{'Color',[0 0 0.9]});
set(hImg,'AlphaData',alphaWeights(2));
uistack(hPlot,'top')

plot([x1(1,:); x2(1,:)],...
     [x1(2,:); x2(2,:)],...
    'Color',[0 0.9 0]);
hold off

if flagSlider
    uicontrol('Style','Slider',...
        'Min',0,'Max',1,'Value',alphaWeights(1),...
        'SliderStep',[0.2 0.2],'Position',[30 10 200 20],...
        'Callback',@(hObj,event) sfm_rawDisplayMatch(img1,img2,x1,x2,'alphaWeights',[get(hObj,'Value') 1-get(hObj,'Value')]))
end
