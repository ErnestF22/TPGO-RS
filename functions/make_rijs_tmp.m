function Rijs_tmp = make_rijs_tmp(rotations_pts, correspondences, num_edges, d, N)
%Make Rijs block matrix, without symmetry-related (i.e. if rotation between
% camera 1 and 2 appears, rotation between 2 and 1 is considered as its inverse
% by default) redundancies

    
%     R=repmat(eye(3),1,1,5)
%     Rinit=rot_randn(R,0.1) %0.1 is the variance (sigma)
    
    Rijs_tmp = zeros(d, d, num_edges);

    for corresp_id = 1:num_edges 
        Rijs_i_id = correspondences(corresp_id, 1);
        Rijs_j_id = correspondences(corresp_id, 1);
            
        if Rijs_i_id > N || Rijs_j_id > N
            disp("Error with correspondences initialization in make_rijs_tmp()");
            break
        end

        %OBS: the step below works (and is setup for) in 3D only i.e. d=3
        angle_diffs_vec = rotations_pts(:, Rijs_j_id) - rotations_pts(:,Rijs_i_id);
        rotz_angle = angle_diffs_vec(3); %Note: angles are expressed in degrees
        roty_angle = angle_diffs_vec(2); %Note: angles are expressed in degrees
        rotx_angle = angle_diffs_vec(1); %Note: angles are expressed in degrees
        Rijs_tmp(:, :, corresp_id) = rotz(rotz_angle) * roty(roty_angle) * rotx(rotx_angle);
    end

end