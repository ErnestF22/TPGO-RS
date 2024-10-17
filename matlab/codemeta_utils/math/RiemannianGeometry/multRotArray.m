%function R=multRotArray(R1,R2)
%Multiply arrays of rotations (R=R1*R2). R1 must be of dimension (3x3xN), R2 must be
%of dimension (3x3xM) where either M==N, N==1 or M==1
function R=multRotArray(R1,R2)
N=size(R1,3);
M=size(R2,3);

if(M~=N && N~=1 && M~=1)
    error('Incorrect array dimensions');
end

N1=max(M,N);

R=zeros(size(R1,1),size(R1,2),N1);

if(N==1)
    for(r=1:N1)
        R(:,:,r)=R1*R2(:,:,r);
    end
else
    if(M==1)
        for(r=1:N1)
            R(:,:,r)=R1(:,:,r)*R2;
        end
    else
        for(r=1:N1)
            R(:,:,r)=R1(:,:,r)*R2(:,:,r);
        end
    end
end
