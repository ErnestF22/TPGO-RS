function GQuery=rotdrd_interpolationLinear(t,G,tQuery)
[R,T]=G2RT(G);
RQuery=rot_interpolationLinear(t,R,tQuery);
TQuery=interp1(t,T',tQuery,'linear')';
GQuery=RT2G(RQuery,TQuery);

