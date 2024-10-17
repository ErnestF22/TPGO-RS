%Compute bearing control field for a single bearing measurement
%function dx=bearingField(y,yg,varargin)
function dx=bearingField(y,yg,varargin)
method='geodesic';
funsAngle=@(x) x^2/2;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'geodesic'
            method='geodesic';
            ivarargin=ivarargin+1;
            funsAngle=lower(varargin{ivarargin});
        case 'projection'
            method='projection';
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

Ny=size(y,2);
Nyg=size(yg,2);
d=size(y,1);
dx=zeros(d,Ny,Nyg);
for iy=1:Ny
    for iyg=1:Nyg
        switch method
            case 'geodesic'
                [v,theta]=cnormalize(sphere_log(y(:,iy),yg(:,iyg)));
                dx(:,iy,iyg)=-funsAngle.f(theta)*v;
            case 'projection'
                Py=orthComplementProjector(y(:,iy));
                dx(:,iy,iyg)=-Py*yg(:,iyg);
        end
    end
end
dx=squeeze(dx);
