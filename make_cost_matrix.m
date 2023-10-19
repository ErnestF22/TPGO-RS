% A = [0 -1; 1 1];
% B = [1 2 3; 4 5 6];
% C = [7 8; 9 10; -11 -12];
% D = [5 5 5; 5 5 5; 5 5 5];
%
% blkmat = [A B; C D];
%
% disp (blkmat)

% R1 = quat2rotm(randrot(1));
% t1 = randn(3, 1);
%
% disp (R1)
% disp (t1)
%
% disp(kron(R1, t1))

n = 10; %number of "frames" -> number of potentially available rot/transl measurements = n^2
d = 3; %dimension of translation/rotation (here: 3D)
%OBS. Changing to d!=3 -> need to generate rotation matrices differently
%(here -> quat2rotm is set up for 3D rotations)

taus = GenerateRandomSymmetricMat(n);
kappas = GenerateRandomSymmetricMat(n);
taus(1:n+1:end) = diag(eye(n)); %set diagonal to 1
kappas(1:n+1:end) = diag(eye(n)); %set diagonal to 1

graph_incidence_mat = GenerateRandomSymmetricBooleanMat(n);
graph_incidence_mat(1:n+1:end) = diag(eye(n)); %set diagonal to 1
% disp(graph_incidence_mat);

taus = taus .* graph_incidence_mat;
kappas = kappas .* graph_incidence_mat;

% R = quat2rotm(randrot(1, n)); %Langevin random would be more appropriate
% R = reshape(R, d, d*n);
% disp(R)
% t = randn(d*n, 1); %this should already be following normal distribution
% disp(t)

upper_triang_size_nodiag = 0.5 * n * (n-1);
upper_triang_size = 0.5 * n * (n+1);


Rnois_nosym = randrot_manopt(d, upper_triang_size_nodiag); %NOTE: NOT all of this will be transformed into rotation matrix
Rnois = eye(d*n);
ii=1;
jj=2;
max_to_reach = n;
% fprintf("size (n-1)n/2 %g\n", size(Rnois_nosym, 1));
for totid = 1:size(Rnois_nosym, 3)
%     Rnois((i-1)*d+1:(i-1)*d+d, (j-1)*d+1:(j-1)*d+d) = quat2rotm(Rnois_nosym(totid));
%     Rnois((j-1)*d+1:(j-1)*d+d, (i-1)*d+1:(i-1)*d+d) = inv(quat2rotm(Rnois_nosym(totid))); % ' instead of inv() would do the same
    Rnois((ii-1)*d+1:(ii-1)*d+d, (jj-1)*d+1:(jj-1)*d+d) = Rnois_nosym(:,:,totid);
    Rnois((jj-1)*d+1:(jj-1)*d+d, (ii-1)*d+1:(ii-1)*d+d) = inv(Rnois_nosym(:,:,totid)); % ' instead of inv() would do the same
    jj = jj+1;
    if (jj>max_to_reach)
        ii = ii+1;
        jj = ii+1;
    end
end
% disp(Rnois)

tnois_part = -5 + (5+5)*rand(3 * 0.5 * n * (n-1), 1); %translation elements are randomized in [-1,1)
tnois = zeros(n*d, n);
% disp(tnois_part);
ii=1;
jj=2;
max_to_reach = n;
for totid = 1:d:size(tnois_part, 1)
    if (ii>=n)
        break
    end
    tnois((ii-1)*d+1:(ii-1)*d+d, jj) = tnois_part(totid:totid+d-1);
    tnois((jj-1)*d+1:(jj-1)*d+d, ii) = -tnois_part(totid:totid+d-1);
    jj = jj+1;
    if (jj>max_to_reach)
        ii = ii+1;
        jj = ii+1;
    end
end
% disp(tnois);

% check if first R matrix has det = 1
% test = [-0.2251   -0.3458   -0.9109; -0.9741    0.0577    0.2188; -0.0231    0.9365   -0.3498];
% disp (det(test));

%%% now forming Q

%L(W^tau)
W_tau_graph = graph(taus);
W_tau_graph_simpl = simplify(W_tau_graph); %probably just puts diag elements to 0
L_Wtau = laplacian (W_tau_graph_simpl); %is this laplacian the correct one?

%L(W^rho)
W_rho_graph = graph(kappas);
W_rho_graph_simpl = simplify(W_rho_graph); %probably just puts diag elements to 0
L_Wrho = laplacian (W_rho_graph_simpl); %is this laplacian the correct one?
% disp(full(L_Wrho)) %!!

%L(Gnois(rho))
L_Gnois_rho = zeros(n*d, n*d);
full_L_Wtau = full(L_Wtau);
for ii = 1:n
    for jj = 1:n
        if ii == jj
            L_Gnois_rho((ii-1)*d+1:(ii-1)*d + d, (jj-1)*d+1:(jj-1)*d + d) = full_L_Wtau(ii, jj) * eye(d);
            continue
        end
        if graph_incidence_mat(ii,jj) ~= 0
            L_Gnois_rho((ii-1)*d+1:(ii-1)*d + d, (jj-1)*d+1:(jj-1)*d + d) = full_L_Wtau(ii, jj) * Rnois((ii-1)*d+1:(ii-1)*d + d, (jj-1)*d+1:(jj-1)*d + d);
            %else
            %matrix already initialized to 0 -> OK
        end
    end
end


%Vnois
Vnois = zeros(n, n*d);
for ii = 1:n
    for jj = 1:n
        if ii == jj
            Vnois(ii, (jj-1)*d+1:(jj-1)*d + d) = sum(taus(ii) .* tnois((jj-1)*d+1:(jj-1)*d + d, :)'); %note: j = i
            continue
        end
        if graph_incidence_mat(jj,ii) ~= 0
            Vnois(ii, (jj-1)*d+1:(jj-1)*d + d) = taus(jj, ii) * tnois((jj-1)*d+1:(jj-1)*d + d, ii); %transpose of tnois() is superfluous in MATLAB
            %else
            %matrix already initialized to 0 -> OK
        end
    end
end


%Sigmanois
Sigmanois = zeros(n*d, n*d);
for diag_id = 1:n %diag_id = i
    sum_tmp = zeros(d);
    for jj = 1:n
        sum_tmp = sum_tmp + (taus(diag_id, jj) * tnois((diag_id-1)*d+1:(diag_id-1)*d + d, jj) * tnois((diag_id-1)*d+1:(diag_id-1)*d + d, jj)');
    end
    Sigmanois(diag_id:diag_id+d-1, diag_id:diag_id+d-1) = sum_tmp;
end


%At last, compute Qnois (and Qnois_tau)
Qnois = L_Gnois_rho + Sigmanois - Vnois' * pinv(full_L_Wtau) * Vnois;
Qnois_tau = Sigmanois - Vnois' * pinv(full_L_Wtau) * Vnois;

disp(Qnois);


