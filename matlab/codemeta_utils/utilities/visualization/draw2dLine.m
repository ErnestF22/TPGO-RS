function draw2dLine(bOrth,varargin)
style='b';
side=1;
plotOptions={};
o=[0;0];

if size(bOrth==3)
    o=homogeneousProjectPoints(o,bOrth);
    bOrth=cnormalize(bOrth(1:2));
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'style'
            ivarargin=ivarargin+1;
            style=varargin{ivarargin};
        case 'side'
            ivarargin=ivarargin+1;
            side=varargin{ivarargin};
        case 'plotoptions'
            ivarargin=ivarargin+1;
            plotOptions=varargin{ivarargin};
        case 'center'
            ivarargin=ivarargin+1;
            o=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

b=computeOrthogonalComplement(bOrth);

c1=o+side*b;
c2=o-side*b;
c=[c1 c2];
plot(c(1,:),c(2,:),style,plotOptions{:})
