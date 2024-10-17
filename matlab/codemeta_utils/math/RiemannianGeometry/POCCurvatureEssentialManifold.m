function POCCurvatureEssentialManifold
I=eye(6);
e3=[0;0;1];
vV=cnormalize([e3;e3]);
bH=[I(:,[1 2 4 5]) [e3;-e3]];

iIdxH=reshape(repmat(1:5,5,1),1,[]);
jIdxH=reshape(repmat(1:5,5,1)',1,[]);

flagNotRepeat=iIdxH>jIdxH;
iIdxH=iIdxH(flagNotRepeat);
jIdxH=jIdxH(flagNotRepeat);

bHBrackets=dblvee(dblbracket(dblhat(bH(:,iIdxH)),dblhat(bH(:,jIdxH))));

disp([(vV'*bHBrackets).^2; iIdxH; jIdxH])



function C=dblbracket(A,B)
C=[ bracket(A(1:3,:,:),B(1:3,:,:));
    bracket(A(4:6,:,:),B(4:6,:,:))];

function C=bracket(A,B)
N=size(A,3);

C=zeros(size(A));

for iA=1:N
    C(:,:,iA)=A(:,:,iA)*B(:,:,iA)-B(:,:,iA)*A(:,:,iA);
end

function vhat=dblhat(v)
N=size(v,2);

vhat=zeros(6,3,N);
for iv=1:N
    vi=v(:,iv);
    vhat(:,:,iv)=[hat(vi(1:3));hat(vi(4:6))];
end

function v=dblvee(vhat)
N=size(vhat,3);

v=zeros(6,N);
for iv=1:N
    vhati=vhat(:,:,iv);
    v(:,iv)=[vee(vhati(1:3,:));vee(vhati(4:6,:))];
end
