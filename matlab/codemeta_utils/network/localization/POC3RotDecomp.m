function POC3RotDecomp

switch 2
    case 1
        %simple test with R2 of the form Rz
        na=sphere_randn(eye(3,1));
        nc=sphere_randn(eye(3,1));

        theta=2*pi*rand;

        Rtheta=blkdiag(1,[cos(theta) -sin(theta); sin(theta) cos(theta)]);
        a=na(2)*nc(2)+na(3)*nc(3);
        b=na(3)*nc(2)-na(2)*nc(3);
        c1=na(1)*nc(1);
        c=na'*Rtheta*nc-c1;
        f=@(cs,sn) a*cos(theta)+b*sin(theta);
        disp(c-f(cos(theta),sin(theta)))

        csEst=(a*c+b*sqrt(a^2+b^2-c^2))/(a^2+b^2);
        snEst=(b*c-a*sqrt(a^2+b^2-c^2))/(a^2+b^2);

        disp(c-f(csEst,snEst))
        disp(csEst^2+snEst^2)
    case 2
        v1=sphere_randn(eye(3,1));
        v2=sphere_randn(eye(3,1));
        v3=sphere_randn(eye(3,1));
        theta1=2*pi*rand;
        theta2=2*pi*rand;
        theta3=2*pi*rand;
        
        R1=rot(v1*theta1);
        R2=rot(v2*theta2);
        R3=rot(v3*theta3);
        R=R1*R2*R3;
        
%         disp([rot(v2*theta2) eye(3)+hat(v2)^2-cos(theta2)*hat(v2)^2+sin(theta2)*hat(v2)]);
%         disp([v1'*R*v3-v1'*v3-v1'*hat(v2)^2*v3 -v1'*hat(v2)^2*v3*cos(theta2)+v1'*hat(v2)*v3*sin(theta2)]);
        a=-v1'*hat(v2)^2*v3;
        b=v1'*hat(v2)*v3;
        c=v1'*R*v3-v1'*v3-v1'*hat(v2)^2*v3;
%        disp([c a*cos(theta2)+b*sin(theta2)]);
        
        csEsta=(a*c+b*sqrt(a^2+b^2-c^2))/(a^2+b^2);
        snEsta=(b*c-a*sqrt(a^2+b^2-c^2))/(a^2+b^2);
        csEstb=(a*c-b*sqrt(a^2+b^2-c^2))/(a^2+b^2);
        snEstb=(b*c+a*sqrt(a^2+b^2-c^2))/(a^2+b^2);

        theta2a=atan2(snEsta,csEsta);
        theta2b=atan2(snEstb,csEstb);
        
        disp('Check that equation for theta2 is valid')
        disp([a*csEsta+b*snEsta-c a*csEstb+b*snEstb-c]);
        disp('Check estimated sine and cosine have correct norm')
        disp([csEsta^2+snEsta^2 csEstb^2+snEstb^2])
        disp('Compare ground truth and recovered theta2')
        disp(modAngle([theta2 theta2a theta2b]))
        
        R2a=rot(v2*theta2a);
        R2b=rot(v2*theta2b);
        
%         disp([R*v3 (eye(3)+hat(v1)^2)*R2*v3+(-cos(theta1)*hat(v1)^2+sin(theta1)*hat(v1))*R2*v3])
%         disp([(R-R2-hat(v1)^2*R2)*v3-[-hat(v1)^2*R2*v3 hat(v1)*R2*v3]*[cos(theta1);sin(theta1)]])
        A1=@(R2) [-hat(v1)^2*R2*v3 hat(v1)*R2*v3];
        b1=@(R2) (R-R2-hat(v1)^2*R2)*v3;
        disp('Check that equation for theta1 is valid')
        disp(max(abs(A1(R2)*[cos(theta1);sin(theta1)]-b1(R2))))
        cssnEsta=A1(R2a)\b1(R2a);
        cssnEstb=A1(R2b)\b1(R2b);
        theta1a=atan2(cssnEsta(2),cssnEsta(1));
        theta1b=atan2(cssnEstb(2),cssnEstb(1));
        disp('Compare ground truth and recovered theta1')
        disp(modAngle([theta1 theta1a theta1b]))

        R1a=rot(v1*theta1a);
        R1b=rot(v1*theta1b);

        A3=@(R2) [-hat(v3)^2*R2'*v1 -hat(v3)*R2'*v1];
        b3=@(R2) (R'-R2'-hat(v3)^2*R2')*v1;
        disp('Check that equation for theta3 is valid')
        disp(max(abs(A3(R2)*[cos(theta3);sin(theta3)]-b3(R2))))
        cssnEsta=A3(R2a)\b3(R2a);
        cssnEstb=A3(R2b)\b3(R2b);
        theta3a=atan2(cssnEsta(2),cssnEsta(1));
        theta3b=atan2(cssnEstb(2),cssnEstb(1));
        disp('Compare ground truth and recovered theta3')
        disp(modAngle([theta3 theta3a theta3b]))

        R3a=rot(v3*theta3a);
        R3b=rot(v3*theta3b);

        disp('Check final solutions')
        disp([R R1a*R2a*R3a;
            R R1b*R2b*R3b])
end