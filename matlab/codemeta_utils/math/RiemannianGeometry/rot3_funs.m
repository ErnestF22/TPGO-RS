function lie_funs=rot3_funs()
lie_funs.dim=3;

lie_funs.exp=@(R,A) expPrivate(R,A);
lie_funs.log=@(R1,R2) logPrivate(R1,R2);

lie_funs.dist=@rot_distn;
lie_funs.mean=@rot_mean;
lie_funs.metric=@rot_metric;

lie_funs.eye=@rot_eye;

lie_funs.parallel=@rot_parallel;

lie_funs.tangentBasis=@rot_tangentBasis;
lie_funs.tangentProj=@rot_tangentProj;
lie_funs.eucl2RiemGrad=@rot_eucl2RiemGrad;
lie_funs.vee=@rot_vee;

lie_funs.randTangentNormVector=@rot_randTangentNormVector;

lie_funs.anchors=@rot_anchors;

lie_funs.retractions=@rot_retractions;

function R2=expPrivate(R,A)
R2=zeros(size(R));
for iR=1:size(A,3)
    R2(:,:,iR)=R*rot(vee(R'*A(:,:,iR)));
end

function A12=logPrivate(R1,R2)
A12=zeros(size(R1));
for iR=1:size(R2,3)
    A12(:,:,iR)=R1*hat(logrot(R1'*R2(:,:,iR)));
end
