function [st,dst,s0,ds0,vVec,ddst,dvVec]=qcrb_randGeodFun(A,varargin)
    [Qct,dQct,Qc0,dQc0,vcVec,ddQct,dvcVec]=rot_randGeodFun(A.Qc,varargin);
    [Rbt,dRbt,Rb0,dRb0,vbVec,ddRbt,dvbVec]=rot_randGeodFun(A.Rb,varargin);
    st=funStruct(Qct,Rbt);
    dst=funStruct(dQct,dRbt);
    s0=valStruct(Qc0,Rb0);
    ds0=valStruct(dQc0,dRb0);
    vVec=valStruct(vcVec,vbVec);
    ddst=funcStruct(ddQct,ddRbt);
    dvVec=valStruct(dvcVec,dvbVec);
end

function f=funcStruct(fQc,fRb)
f=@(t) struct('Qc',fQc(t),'Rb',fRb(t));
end

function v=valStruct(Qc,Rb)
v=struct('Qc',Qc,'Rb',Rb);
