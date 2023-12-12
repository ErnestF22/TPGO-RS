function [a,b]=hingeSparseBoost(x,y,w,v,varargin)
C=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'c'
            ivarargin=ivarargin+1;
            C=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NPoints=size(x,2);
NLearners=size(w,2);
r=hingeResponse(x,w,v,'appendOpposite');
u=ones(NPoints,1);

cvx_begin
    variables a(2*NLearners,1) b z(NPoints,1)
    minimize (norm(w,1)+C*sum(z))
    subject to
        y.*(r*a+b)-u+z>=0
        z>=0
cvx_end
