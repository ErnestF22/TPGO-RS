function test_recovery_recap


N = 6;
load('Qbnn_data/testdata.mat','edges')
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')


load('Qbnn_data/testdata.mat','R')
load('Qbnn_data/testdata.mat','T')
load('Qbnn_data/testdata.mat','Tijs')

nrs = size(R, 1);
p = nrs;

d = size(R, 2);

low_deg = 2;

node_degrees = [2,2,4,3,3,4];
nodes_high_deg = node_degrees == 3 | node_degrees == 4;
nodes_low_deg = ~nodes_high_deg;


T_edges = make_T_edges(T, edges);

RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)), T_edges];
Qx = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
RT_stacked_high_deg_poc = Qx * RT_stacked_high_deg;

R_tilde = multiprod(repmat(Qx, 1, 1, N), R);


disp('RT_stacked_high_deg_poc')
disp(RT_stacked_high_deg_poc)

disp("max(RT_stacked_high_deg_poc(d+1:end, :), [], ""all"")")
disp(max(RT_stacked_high_deg_poc(d+1:end, :), [], "all"))

disp('matStackH(R_recovered)')
disp(matStackH(R_tilde))

%%
% Our goal is to find Rbi_hat that compensates Qbi s.t.
% P_last*Qbi*Ri_manoptP_last*blkdiag(eye, Rbi_hat)*Ri_manopt = 0
R_tilde_low_deg = R_tilde(:,:,nodes_low_deg);

P_last = make_p_last(p,d);

%% now focus on node ii=1
iis = [1 2];

for ii = iis
    disp("Running full recovery procedure on node:");
    disp(ii)
    Ri_tilde = R_tilde_low_deg(:,:,ii);

    %% preliminary checks

    %computing Tij12
    Tij1j2 = zeros(d,low_deg);
    num_js_found = 0;
    for e = 1:size(edges)
        e_i = edges(e,1);
        %         e_j = edges(e,2);
        if (e_i == ii)
            num_js_found = num_js_found + 1;
            Tij1j2(:,num_js_found) = Tijs(:,e);
        end
    end

    %checking eq. (61)
    check_61 = Ri_tilde * Tij1j2;
    tail_check_61 = check_61(end,:); %TODO: make this smarter using tail()
    disp("max(tail(check_61, p-d), [], ""all"")");
    disp(max(tail_check_61, [], "all"));

    %checking eq. (62)
    check_62 = P_last * Ri_tilde * Tij1j2;
    disp("max(check_62, [], ""all"")");
    disp(max(check_62, [], "all"));

    %% 3-step procedure

    % step 1
    Ri_tilde_2 = Ri_tilde(size(Ri_tilde, 1) - (p - 2) + 1:end, :);
    disp('svd(Ri_tilde_2)')
    disp(svd(Ri_tilde_2))
    % step 2
    R_hat_b_i_last = orthCompleteBasis(Ri_tilde_2);
    % step 3
    Rb_i = orthCompleteBasis(R_hat_b_i_last);

    %% final checks
    Qb_i = blkdiag(eye(p-node_deg), Rb_i); %node_deg = 2

    %checking eq. (63)
    check_63 = Qb_i * Ri_tilde * Tij1j2;
    tail_check_63 = check_63(end,:); %TODO: make this smarter using tail()
    disp("max(tail(check_63, p-d), [], ""all"")");
    disp(max(tail_check_63, [], "all"));

    %checking eq. (66)
    check_66 = (Qx' * Ri_tilde) * Tij1j2;
    tail_check_66 = check_66(end,:); %TODO: make this smarter using tail()
    disp("max(tail(check_66, p-d), [], ""all"")");
    disp(max(tail_check_66, [], "all"));

    %checking eq. (67)
    check_67 = P_last * blkdiag(eye(p-d), R_hat_b_i_last) * Ri_tilde;
    disp("max(check_67, [], ""all"")");
    disp(max(check_67, [], "all"));

    %checking eq. (68)
    %     check_68 = R_hat_b_i_last * Ri_tilde_2;
    %     disp("max(check_68, [], ""all"")");
    %     disp(max(check_68, [], "all"));


    % "final test"
    Ri_stief = blkdiag(eye(p-d), Rb_i) * Ri_tilde;
    disp("Ri_stief");
    disp(Ri_stief);
    for e = 1:size(edges)
        e_i = edges(e,1);
        e_j = edges(e,2);
        if (e_i == ii)
            disp("e_i, e_j")
            disp([e_i, e_j])
            m_check = Ri_stief * Tijs(:,e) - T(:,e_j) + T(:,e_i);
            disp("max(m_check, [], ""all"")");
            disp(max(m_check, [], "all"));
        end
    end



end

% check if recovery has actually been performed correctly
% ...




end %file function