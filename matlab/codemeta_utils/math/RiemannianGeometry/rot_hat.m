function A=rot_hat(R,a)
NR=size(R,3);
Na=size(a,2);
d=size(R,1);
if Na>1
    A=zeros(d,d,Na);
    if NR>1
        if Na~=NR
            error('If R contains more than one rotation, A should contain the same number of tangent vectors')
        end
        for n=1:Na
            A(:,:,n)=rot_hat(R(:,:,n),a(:,n));
        end
    else
        for n=1:Na
            A(:,:,n)=rot_hat(R,a(:,n));
        end
    end
else
    A=zeros(d,d);
    cnt=1;
    %the sign and order give consistency with the hat() and vee() functions for SO(3)
    for c=d-1:-1:1
        for r=d:-1:c+1
            A(r,c)=(-1)^(r+c+1)*a(cnt);     
            cnt=cnt+1;
        end
    end
    A=A-A';
    A=R*A;
end
