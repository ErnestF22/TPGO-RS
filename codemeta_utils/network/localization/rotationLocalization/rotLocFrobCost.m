function [c,gradc]=rotLocFrobCost(E,R,RRel,varargin)
[c,gradc]=rotLocBaseCost(E,R,RRel,@rotLocFrobCostPair,varargin{:});
