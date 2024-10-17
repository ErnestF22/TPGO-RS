function essential_flipAmbiguity_test
switch 2
    case 1
        Q=essential_randn();
        E=zeros(3,3,4);
        for k=1:4
            E(:,:,k)=essential_getE(essential_flipAmbiguity(Q,k));
        end
        E=reshape(E,3,[]);
        disp([E; E-kron([1 1 -1 -1],E(1:3,1:3))])
    case 2
        Q=zeros(6,3);
        Q(1:3,:)=rot(3*pi/8*[0;1;0]);
        Q(4:6,:)=rot(5*pi/8*[0;1;0]);
        for k=1:4
            QFlip=essential_flipAmbiguity(Q,k);
            hold on
            plot3(2,0,1.5,'*')
            hold off
            
            subplot(2,2,k)
            essential_disp(QFlip)
            axis square
            axis equal
        end
    case 3
        
        
end
