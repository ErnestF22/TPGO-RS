function M = get_marix_M(E,F,nV,tiji,tiji_mtrx)




M =zeros(3*nV,3*nV);

for ie=1:size(E,1)
    

    i   =  E(ie,1);
    j  =  E(ie,2);
    
    idxr = 3*i-2:3*i;
    idxc = 3*j-2:3*j;
    
    M(idxr,idxc) =  tiji(:,1,ie)*tiji(:,2,ie)';
    
        for kk =1:numel(F{ie})
            
            k = F{ie}(kk);
            nijki = cross(  tiji_mtrx(3*i-2:3*i,j),  tiji_mtrx(3*i-2:3*i,k));
            nijki = nijki / norm(nijki);
            
            njikj = cross(  tiji_mtrx(3*j-2:3*j,i),  tiji_mtrx(3*j-2:3*j,k));
            njikj = njikj / norm(njikj);
            
            M(idxr,idxc) = M(idxr,idxc) + nijki*njikj';
        end
   
    
end

M = M + M';
 

 


end

