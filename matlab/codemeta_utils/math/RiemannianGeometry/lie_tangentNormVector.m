function v=lie_tangentNormVector(lf,y,v)
v=v/sqrt(lf.metric(y,v,v));
