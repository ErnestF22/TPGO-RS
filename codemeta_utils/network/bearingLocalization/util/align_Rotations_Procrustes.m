function [Yhat_al_3d,err_deg,err_rad] = align_Rotations_Procrustes(Rgt,Rb)


nV = size(Rgt,3);
Ygt= zeros(3,3*nV);
Yhat= zeros(3,3*nV);
Yhat_al_3d = zeros(size(Rb));
for i=1:nV,Ygt(:,3*i-2:3*i) = Rgt(:,:,i); end 
for i=1:nV,Yhat(:,3*i-2:3*i) = Rb(:,:,i); end
[U,~,V]  = svd(Yhat*Ygt');
Ral = V*diag([1 1 det(V*U')])*U';
Yhat_al = Ral*Yhat;
for i=1:nV,Yhat_al_3d(:,:,i) = Yhat_al(:,3*i-2:3*i);  end



% now measure the errors in degrees
err_rad = zeros(nV,1);
for i=1:nV
    err_rad(i) = acos(.5*trace(Rgt(:,:,i)'*Yhat_al_3d(:,:,i)) -0.5);
end
err_deg = err_rad*180/pi;
 



end

