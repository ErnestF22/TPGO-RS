function householderRotation_test
D=5;
k=4;
E=eye(D);

v1=randn(D,1);
for mode=1:2
    switch mode
        case 1
            v2=E(:,k);
            H=householderRotation(v1,k);
        case 2
            v2=randn(D,1);
            H=householderRotation(v1,v2);
    end
    display(H)

    disp('cnormalize([v2 H*v1])')
    disp(cnormalize([v2 H*v1]))

    disp('det(H)')
    disp(det(H))

    disp('H''*H-eye(D)')
    disp(H''*H-eye(D))
end
