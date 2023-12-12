function POCPosdefHouseholder
%find a pd matrix S that sends w to v, where v'*w>0
%resetRands(2)

d=3;
v=randn(d,1);

while 1
    w=randn(d,1);
    
    if v'*w>0
        break
    end
end

vNorm=cnormalize(v);

p=vNorm'*w;
wOrth=w-p*vNorm;
[wOrthNorm,s]=cnormalize(wOrth);

R=orthCompleteBasis([vNorm wOrthNorm])';

avw=atan2(s,p);

aAdj=pi/4-avw/2;
RAdj=rot2(aAdj);
RAdj=blkdiag(RAdj,eye(d-2));

M=RAdj*R*[v w];
%disp(RAdj*R*[vNorm cnormalize(w)])

L=blkdiag(diag(M(1:2,1)./M(1:2,2)),eye(d-2));
disp(diag(L))

S=R'*(RAdj'*L*RAdj)*R;
disp(v-S*w)
disp(S-S')


function R=rot2(a)
R=[cos(a) -sin(a); sin(a) cos(a)];
