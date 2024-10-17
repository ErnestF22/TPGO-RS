function [REst,TEst]=poseEstimationRefineFromRT(RInit,TInit,X,x,varargin)
GEst=poseEstimationRefineFromG(RT2G(RInit,TInit),X,x,varargin{:});
REst=G2R(GEst);
TEst=G2T(GEst);
