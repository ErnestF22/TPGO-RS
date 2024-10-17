function [A,b,R_gf] = make_A_b(R_gf, T_globalframe_nois, Tijs_vec, edges, params)
%MAKE_A_B function returns A,b matrices used in step 2 of Manopt-related
%translation estimation procedure
%and covert R to "3D" mode in case it is passed as stacked

if ~exist('params','var') || isempty(params)
    N = size(R_gf, 3);
    d = size(R_gf, 2);
    num_edges = size(edges,1);
else
    N = params.N;
    d = params.d;
    d_aff = params.d_aff;
    global_camera_id = params.global_camera_id;
    num_tests_per_sigma = params.num_tests_per_sigma;
    transf_end_thresh = params.transf_end_thresh;
    max_icp_iterations = params.max_icp_iterations;
    num_edges_full = params.num_edges_full;
    num_edges = size(edges,1);
    procrustes_mode = params.procrustes_mode;
    initguess_is_available = params.initguess_is_available;
end

idxEdges=reshape(1:d*num_edges,d,num_edges);
idxNodes=reshape(1:d*N,d,N);

if (size(R_gf, 3)>1)
    A=zeros(d*num_edges,d*N);
    b=zeros(d*num_edges, 1);
%     b = Tijs_vec(:);    
    for edge_id = 1:num_edges
        ii = edges(edge_id, 1);
        jj = edges(edge_id, 2);
        %             A(col_id*d-(d-1):col_id*d, (ii*d)-(d-1):ii*d) = R(:,:,ii)';
        %             A(col_id*d-(d-1):col_id*d, (jj*d)-(d-1):jj*d) = -R(:,:,ii)';
        A(idxEdges(:,edge_id),idxNodes(:,ii)) = R_gf(:,:,ii)';
        A(idxEdges(:,edge_id),idxNodes(:,jj)) = -R_gf(:,:,ii)';
        b(idxEdges(:,edge_id)) = Tijs_vec(:, edge_id);
    end
else
    R2 = zeros(d,d,N);
    A=zeros(d*num_edges,d*N);
    b=zeros(d*num_edges, 1);
%     b = Tijs_vec(:);
    for edge_id = 1:num_edges
        ii = edges(edge_id, 1);
        jj = edges(edge_id, 2);
        %             A(col_id*d-(d-1):col_id*d, (ii*d)-(d-1):ii*d) = R(:,ii*d-(d-1):ii*d)';
        %             A(col_id*d-(d-1):col_id*d, (jj*d)-(d-1):jj*d) = -R(:,ii*d-(d-1):ii*d)';
        A(idxEdges(:,edge_id),idxNodes(:,ii)) = R_gf(ii*d-(d-1):ii*d, :)';
        A(idxEdges(:,edge_id),idxNodes(:,jj)) = -R_gf(ii*d-(d-1):ii*d, :)';
        b(idxEdges(:,edge_id)) = Tijs_vec(:, edge_id);
    end
    for ii = 1:N
        R2(:,:,ii) = R_gf(ii*d-(d-1):ii*d, :); %covert R to "3D mode"
    end
    R_gf = R2;
end

end

