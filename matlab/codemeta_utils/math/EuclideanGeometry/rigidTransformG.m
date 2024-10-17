function [XTransformed,JXTransformed,HXTransformed]=rigidTransformG(G,X,varargin)
R=G2R(G);
T=G2T(G);
switch nargout
    case {0,1}
        XTransformed=rigidTransform(R,T,X,varargin{:});
    case 2
        [XTransformed,JXTransformed]=rigidTransform(R,T,X,varargin{:});
    case 3
        [XTransformed,JXTransformed,HXTransformed]=rigidTransform(R,T,X,varargin{:});
end