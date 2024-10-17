function [tijgt,tijigt,tiji,tiji_mtrx] = generate_baerings(Rgt,Tgt,E,nV,sigma)


nE  = size(E,1);


manifold = spherefactory(3);


tijgt = zeros(3,2,nE) ; % [tij tji] in the global reference frame
tijigt = zeros(3,2,nE) ; % [tij tji] in the  local reference frame
tiji = zeros(3,2,nE); % noisy local baerings
tiji_mtrx = zeros(3*nV,nV); 

for ie=1:nE
   
    i  =  E(ie,1);
    j  =  E(ie,2);
    
    tijgt(:,1,ie) = (Tgt(:,j) - Tgt(:,i))/norm(Tgt(:,j) - Tgt(:,i));
    tijgt(:,2,ie) = (Tgt(:,i) - Tgt(:,j))/norm(Tgt(:,j) - Tgt(:,i));
    
    tijigt(:,1,ie) = Rgt(:,:,i)'*tijgt(:,1,ie);
    tijigt(:,2,ie) = Rgt(:,:,j)'*tijgt(:,2,ie);
       
    tiji(:,1,ie) = manifold.exp(  tijigt(:,1,ie), sigma*randn*manifold.randvec(tijigt(:,1,ie))   );
    tiji(:,2,ie) = manifold.exp(  tijigt(:,2,ie), sigma*randn*manifold.randvec(tijigt(:,2,ie))   );
    
    tiji_mtrx(3*i-2:3*i,j) =   tiji(:,1,ie);
    tiji_mtrx(3*j-2:3*j,i) =   tiji(:,2,ie);
    
end


end