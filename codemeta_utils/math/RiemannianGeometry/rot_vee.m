function a=rot_vee(R,A)
if(size(A,3)>1)
    if(size(R,3)>1)
        if(size(R,3)~=size(A,3))
            error('If R contains more than one rotation, A should contain the same number of tangent vectors')
        end
        for n=1:size(R,3)
            a(:,n)=rot_vee(R(:,:,n),A(:,:,n));
        end
    else
        for n=1:size(A,3)
            a(:,n)=rot_vee(R,A(:,:,n));
        end
    end
else
    A=R'*A;
    n=size(A,1);
    a=zeros(n*(n-1)/2,1);
    cnt=1;
    %the sign and order give consistency with the hat() and vee() functions for SO(3)
    for c=n-1:-1:1
        for r=n:-1:c+1
            a(cnt)=(-1)^(r+c+1)*A(r,c);     
            cnt=cnt+1;
        end
    end
end
