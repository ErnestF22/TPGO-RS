function [eulerAngles]=rotationMatrix2EulerAngles(rotationMatrix)

if(rotationMatrix(3,1)~=1&&rotationMatrix(3,1)~=-1)
    
    theta1=-asin(rotationMatrix(3,1));
    theta2=pi-theta1;
    phi1=atan2(rotationMatrix(3,2)/cos(theta1),rotationMatrix(3,3)/cos(theta1));
    phi2=atan2(rotationMatrix(3,2)/cos(theta2),rotationMatrix(3,3)/cos(theta2));
    psi1=atan2(rotationMatrix(2,1)/cos(theta1),rotationMatrix(1,1)/cos(theta1));
    psi2=atan2(rotationMatrix(2,1)/cos(theta2),rotationMatrix(1,1)/cos(theta2));
    
    eulerAngles=[phi1 phi2;theta1 theta2;psi1 psi2];
     
else
    
    psi=0;
    
    if(rotationMatrix(3,1)==-1)
        theta=pi/2;
        phi=psi+atan2(rotationMatrix(1,2),rotationMatrix(1,3));
    else
        theta=-pi/2;
        phi=-psi+atan2(-rotationMatrix(1,2),rotationMatrix(1,3));
    end
    
    eulerAngles=[phi phi;theta theta;psi psi];

end