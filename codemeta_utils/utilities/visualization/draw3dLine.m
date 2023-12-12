function draw3dLine(bOrth,varargin)
style='b';
side=1;
plotOptions={};
o=[0;0;0];

if size(bOrth,1)==4
    bOrth=planeNormalize(bOrth);
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

if size(bOrth,1)==4
    o=homogeneousProjectPoints(o,bOrth);
    bOrth=bOrth(1:3,:);
end

[U,~,~]=svd(bOrth);
b=U(:,3);

c1=o+side/2*b;
c2=o-side/2*b;
c=[c1 c2];
plot3(c(1,:),c(2,:),c(3,:),style,plotOptions{:})
