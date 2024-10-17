function draw3dcameraFromG(G,varargin)
[R,T]=G2RT(G);
draw3dcameraFromRT(R,T,varargin{:});