function [v,Q2r]=essential_log(Q1,Q2,varargin)
NQ=size(Q2,3);
Q2r=essential_closestRepresentative(Q1,Q2,varargin{:});
v=zeros(6,3,NQ);
for iQ=1:NQ
    v(1:3,:,iQ)=rot_log(Q1(1:3,:),Q2r(1:3,:,iQ));
    v(4:6,:,iQ)=rot_log(Q1(4:6,:),Q2r(4:6,:,iQ));
end

