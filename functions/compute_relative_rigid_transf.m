function transf_ij = compute_relative_rigid_transf(transf_i, transf_j)
%Takes as input two affine transformations i, j and returns the
%relative transform that links i to j
    
    %     transf_ij = zeros(size(transf_i)); % =
    %     zeros(size(transf_j))
    % it is overwritten anyways
    
    disp("NOTE: compare_relative_aff_transf() works only with 3D rigid transforms!");

    rigid_i = rigidtform3d(transf_i);
    rigid_j = rigidtform3d(transf_j);

    rot_ij = rigid_j.R * (rigid_i.R)';

    transl_ij = rigid_j.Translation - rigid_i.Translation;

    transf_ij = rigidtform3d(rot_ij,transl_ij).A;

end %function