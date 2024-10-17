function GQuery=rotdrd_interpolationProjected(t,G,tQuery,varargin)
[R,T]=G2RT(G);
RQuery=rot_interpolationProjected(t,R,tQuery,varargin{:});
TQuery=real_interp(t,T,tQuery,varargin{:});
GQuery=RT2G(RQuery,TQuery);

