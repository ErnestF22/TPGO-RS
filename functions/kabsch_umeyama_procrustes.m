function transf_curr = som_procrustes(T_globalframe_nois, Tijs_vec, edges, params)
% SOM_SIMPLE_COORD_DESC_PROCRUSTES Computes transformation through 
% Procrustes pipeline
% Inputs can obviously be noisy

%copy params
N = params.N;
d = params.d;
d_aff = params.d_aff;
global_camera_id = params.global_camera_id;
num_tests_per_sigma = params.num_tests_per_sigma;
transf_end_thresh = params.transf_end_thresh;
max_icp_iterations = params.max_icp_iterations;
num_edges_full = params.num_edges_full;
num_edges = params.num_edges;
procrustes_mode = params.procrustes_mode;
initguess_is_available = params.initguess_is_available;


% dummy Kabsch
% a = [1 1; 2 2];
% b = [2 2; 3 3];
% [U,r,lrms] = Kabsch(a,b)

Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

%initialize initial guesses (if gt available -> use it here)
eye_so_dn = ones(d*N, d);
for ii = 1:d:d*N
    eye_so_dn(ii:ii+d-1, :) = eye(d);
end
eye_so_dn = matUnstack(eye_so_dn);

transl_zero_n = zeros(d*N, 1);

rot_prev =  eye_so_dn;
transl_prev = transl_zero_n;
transf_prev = make_transf(rot_prev, transl_prev);

R = zeros(N*d,d);


%Procrustes - step 1: Kabsch/Umeyama
for ii = 1:N
    X = zeros(d, N-1); %following Umeyama notation
    Y = zeros(d, N-1);
    for col_id = 1:size(edges, 2)
        if (edges(1, col_id) == ii)
            jj = edges(2, col_id);
        end
        if (edges(2, col_id) == ii)
            jj = edges(1, col_id);
        end
        if jj > ii
            xy_mats_id = jj-1;
        else
            xy_mats_id = jj;
        end
        X (:, xy_mats_id) = Tijs_vec(:,col_id);
        %         Y (:, xy_mats_id) = T(:, j) - T(:,i);
        Y (:, xy_mats_id) = Tijs_vec(:,col_id);
    end
    %     disp(X)
    %     disp(Y)
    %             transform = procrustes_kabsch(T(:, j) - T(:,i), Tijs(:,col_id), d);
    procrustes_output_R = eye(d);
    if(procrustes_mode == "kabsch")
        procrustes_output_R = procrustes_kabsch(X,Y,d);
        R((ii-1)*d+1:(ii-1)*d+d, :) = procrustes_output_R;
    elseif (procrustes_mode == "umeyama")
        procrustes_output_R = procrustes_umeyama(X,Y,d).R;
        R((ii-1)*d+1:(ii-1)*d+d, :) = procrustes_output_R;
    else
        disp("Error: Procrustes mode unknown! Using default");
        R = som_stepone_procrustes(T_globalframe_nois, Tijs_vec, Tijs_mat, edges, params);
        break;
    end
end

%COORD DESC - step 2
% T = euclideanfactory(d, N);
% step2.M = T;
% step2.cost = @(x) norm(A * T(:) + b);
% [T, xcost2, info2, options2] = trustregions(step2);

R = matUnstack(R);
[T,A,b] = som_steptwo_procrustes(R, T_globalframe_nois, Tijs_vec, Tijs_mat, edges, params);
% T = reshape(T, d, []);

transf_curr = make_transf(R, T);

%WHEN ITERATING... UNCOMMENT BELOW
% num_iterations = 0;
% %COORD DESC - step 3: iterate until convergence
% while (norm(transf_prev - transf_curr)>= transf_end_thresh && num_iterations<max_icp_iterations)
%     rot_prev =  R;
%     transl_prev = T;
%     transf_prev = transf_curr;
% 
% 
%     %COORD DESC - step 1
%     R = zeros(N*d,d);
%     for ii = 1:N
%         X = zeros(d, N-1); %following Umeyama notation
%         Y = zeros(d, N-1);
%         for col_id = 1:size(correspondences, 2)
%             if (correspondences(1, col_id) == i)
%                 jj = correspondences(2, col_id);
%     
%             end
%             if (correspondences(2, col_id) == i)
%                 jj = correspondences(1, col_id);
%             end
%             if jj > i
%                 xy_mats_id = jj-1;
%             else
%                 xy_mats_id = jj;
%             end
%             X (:, xy_mats_id) = Tijs_tmp_nois(:,col_id);
%             %         Y (:, xy_mats_id) = T(:, jj) - T(:,ii);
%             Y (:, xy_mats_id) = Tijs_tmp_nois(:,col_id);
%         end
%         %     disp(X)
%         %     disp(Y)
%         %             transform = procrustes_kabsch(T(:, jj) - T(:,ii), Tijs(:,col_id), d);
%         if(procrustes_mode == "kabsch")
%             procrustes_output_R = procrustes_kabsch(X,Y,d);
%             R((ii-1)*d+1:(ii-1)*d+d, :) = procrustes_output_R;
%         elseif (procrustes_mode == "umeyama")
%             procrustes_output_R = procrustes_umeyama(X,Y,d).R;
%             R((ii-1)*d+1:(ii-1)*d+d, :) = procrustes_output_R;
%         else
%             disp("Error: Procrustes mode unknown!");
%             R = som_stepone_procrustes(T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, params);
%             break;
%         end
%         procrustes_output_R = procrustes_umeyama(X,Y,d);
%         R((ii-1)*d+1:(ii-1)*d+d, :) = procrustes_output_R.R;
%     end
%     
% 
%     %COORD DESC - step 2
%     [T,A,b] = som_steptwo_procrustes(R, T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, N, d);
%     % T = reshape(T, d, []);
% 
%     transf_curr = make_transf(R,T);
% end
% 
% 
% % disp(T);
% disp(transf_curr);

end %function
