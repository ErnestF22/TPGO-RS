%function [Rt,vt]=rot_geodFun(R,v)
%Generate functions for the rotations Rt(t) and tangent vectors vt(t) that
%describe the geodesic rot_exp(R,t*v). If R is omitted or empty, generate a
%random R. If v is omitted, generate a random normal geodesic starting from
%R.
%
%See also rot_randGeodFun
function [Rt,dRt,R0,dR0,vVec,ddRt,dvVec]=qcrb_geodFun(R0)

    [Rt.qc,dRt.qc,R0.qc,dR0.qc,vVec.qc,ddRt.qc,dvVec.qc] = rot_geodFun(R0.qc, []);
    [Rt.rb,dRt.rb,R0.rb,dR0.rb,vVec.rb,ddRt.rb,dvVec.rb] = rot_geodFun(R0.rb, []);

end % file function

