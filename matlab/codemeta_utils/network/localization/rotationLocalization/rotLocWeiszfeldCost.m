function [c,gradc]=rotLocWeiszfeldCost(E,R,RRel,varargin)
if nargout==1
    c=rotLocBaseCost(E,R,RRel,@rotLocWeiszfeldCostPair,varargin{:});
else
    [c,gradc]=rotLocBaseCost(E,R,RRel,@rotLocWeiszfeldCostPair,varargin{:});
end
