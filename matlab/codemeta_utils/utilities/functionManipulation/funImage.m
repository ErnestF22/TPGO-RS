%Show function as an image
%function z=funImage(x,y,f,varargin)
%Evaluates the function f([x;y]) on meshgrid(x,y) and shows the result
%using imagesc
function z=funImage(x,y,f,varargin)
method='image';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

lx=length(x);
ly=length(y);

z=zeros(lx,ly);

for ix=1:lx
    for iy=1:ly
        z(ix,iy)=f([x(ix);y(iy)]);
    end
end

switch method
    case 'image'
        imagesc(x,y,z')
        set(gca,'YDir','normal')
    case 'surf'
        surf(x,y,z','EdgeColor','none');
        view(0,90)
    case 'contour'
        contour(x,y,z',100);
    otherwise
        error('Invalid imagefun method');
end
