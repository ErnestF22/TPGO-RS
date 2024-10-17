%Show indicator set of 2-dimensional symbolic constraints
%function funImageConstraints(constraints,x,y,params,varargin)
%Arguments
%   constraints     array of inequality constraints
%   x,y             the two independent variables to plot
%   params          cell array of substitutions of the form
%       {var,val,var,val,...}
%Optional arguments
%   'xgrid',xgrid       x values for the evaluation grid
%   'ygrid',ygrid       y values for the evaluation grid
%   'optsFunImage',opts options for call to funImage
function funImageConstraints(constraints,x,y,params,varargin)
xGrid=linspace(0,1,50);
yGrid=linspace(0,1,50);
optsFunImage={};
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'xgrid'
            ivarargin=ivarargin+1;
            xGrid=varargin{ivarargin};
        case 'ygrid'
            ivarargin=ivarargin+1;
            yGrid=varargin{ivarargin};
        case 'optsfunimage'
            ivarargin=ivarargin+1;
            optsFunImage=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

constraints=subs(constraints,params(1:2:end),params(2:2:end));

f=@(v) all(subs(constraints,{x y},{v(1) v(2)}));
funImage(xGrid,yGrid,f,optsFunImage{:});
