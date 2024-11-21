function [xt,dxt]=somLambda_geodFun(x0,vx0)
[lambdat,dlambdat]=real_geodFun(x0.lambda,vx0.lambda);
[Tt,dTt]=real_geodFun(x0.T,vx0.T);
[Rt,dRt]=rot_geodFun(x0.R,vx0.R);

function xtEval=xtFun(t)
xtEval=struct('lambda',lambdat(t),'T',Tt(t),'R',Rt(t));
end
function dxtEval=dxtFun(t)
dxtEval=struct('lambda',dlambdat(t),'T',dTt(t),'R',dRt(t));
end

xt=@xtFun;
dxt=@dxtFun;

end
