function [st,dst,s0,ds0,vVec,ddst,dvVec]=qcrb_randGeodFun(A,varargin)
    [qct,dqct,qc0,dqc0,vcVec,ddqct,dvcVec]=rot_geodFun(A.qc,varargin);
    [rbt,drbt,rb0,drb0,vbVec,ddrbt,dvbVec]=rot_geodFun(A.rb,varargin);
    st=funcStruct(qct,rbt);
    dst=funcStruct(dqct,drbt);
    s0=valStruct(qc0,rb0);
    ds0=valStruct(dqc0,drb0);
    vVec=valStruct(vcVec,vbVec);
    ddst=funcStruct(ddqct,ddrbt);
    dvVec=valStruct(dvcVec,dvbVec);
end

function f=funcStruct(fqc,frb)
f=@(t) struct('qc',fqc(t),'rb',frb(t));
end

function v=valStruct(qc,rb)
v=struct('qc',qc,'rb',rb);
end
