function transf_ij = compute_relative_aff_transf(transf_i, transf_j)
%Takes as input two affine transformations i, j and returns the
%relative transform that links i to j
    
%     transf_ij = zeros(size(transf_i)); % = zeros(size(transf_j))
    
    rot_ij = G2R(transf_j) * G2R(transf_i)';

    transl_ij = G2T(transf_j) - G2T(transf_i);

    transf_ij = RT2G(rot_ij,transl_ij);

end