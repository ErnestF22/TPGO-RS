function retval = procrustes_kabsch(a,b, dim)
%Compute rotation according to Kabsch's Procrustes formulation
    Tijs_kabsch = a';
    Tijs_true_kabsch = b';
    
    mean_Tijs = mean(Tijs_kabsch); %comes as row vector!
    mean_Tijs_true = mean(Tijs_true_kabsch); %comes as row vector!
    
    Tijs_kabsch = Tijs_kabsch - mean_Tijs;
    Tijs_true_kabsch = Tijs_true_kabsch - mean_Tijs_true;
    
    H = Tijs_kabsch' * Tijs_true_kabsch;
    
    % R = (mpower(H'*H,0.5)) * inv(H); %works only when H is invertible
    
    [U,S,V] = svd(H);
    detsign = sign(det(V*U'));
    S_R = eye(dim);
    S_R(dim,dim) = detsign;
    R = V * S_R * U';
    retval = R;
end