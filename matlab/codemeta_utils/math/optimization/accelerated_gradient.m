% Accelerated Gradient Method
%function accelerated_gradient(gradient, y0, alpha, maxIt)
% Inputs
%   gradient  Function to compute gradient
%   y0        Initial value for x
%   alpha     Step size
%   maxIt     Max number of iterations
% Output
%   xs        Values for all iterates
function xs=accelerated_gradient(gradient, y0, alpha, maxIt)
ys=[y0 zeros(1,maxIt)];
xs=[y0 zeros(1,maxIt)];
for t=1:(maxIt+1)
    y=ys(t);
    x=xs(t);
    g=gradient(y);
    x_plus=y-alpha*g;
    y_plus=x_plus+((t-1)/(t+2))*(x_plus-x);
    ys(t+1)=y_plus;
    xs(t+1)=x_plus;
end
