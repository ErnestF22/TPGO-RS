function H=groundHomographyFromRT(R1,T1,varargin)
H=groundHomographyFromG(RT2G(R1,T1),varargin{:});
