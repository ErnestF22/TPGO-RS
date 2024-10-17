function [Y0, lambda_pim_out, v_pim_out] = rsom_pim_hessian( ...
    R, problem_struct_next, thresh)
%RSOM_PIM_HESSIAN Return a new starting point Y0 with lower cost that R
% This is based on a linesearch towards an eigenvector v_pim_out 
% corresponding to negative eigenvalue lambda_pim_out.
% If the Hessian does not have any negative eigenvalue (i.e., the two PIM
% iterations do not find it), simply return Y0 = R and the maximum
% eigenvalue with an associated eigenvector.

x = cat_zero_rows_3d_array(R);
rhess_fun_han = @(u) rsom_rhess_rot_stiefel(x,u,problem_struct_next);

stiefel_normalize_han = @(x) x./ (norm(x(:)));

u_start = stiefel_randTangentNormVector(x);
% u_start = [0.0163755511583415	0.047206139253622	0.043744183451971	0.0164468525265251	-0.183939866321754	0.0775765615568843	0.081705893500047	-0.0763981114818336	-0.0681940471296265	0.0200066632490414	-0.0693327179989984	-0.0100439354072153	0.0796844674113843	-0.00900207458640181	-0.0254696760262091;
% -0.0517352939056209	0.00255667669961411	0.00102901731330895	0.000121213492870357	0.0647676867712249	0.0576003051088084	0.125060573308446	-0.1344844404486	0.0571221077372798	0.00471263201973388	-0.105468994360747	-0.0232746551769924	0.0257593997615303	-0.00457408825517136	0.0791961808260063;
% -0.0412363031609128	0.01696780405102	0.0145233439917445	-0.00891040373179798	-0.0596906867066129	-0.179455806586354	0.0208563340870208	-0.000148083566842182	-0.195523524716021	-0.122258885953187	-0.00672241060642059	-0.0395823480979852	0.00190089304860734	3.80393074642591e-05	-0.0138893627510208;
% -0.19090592185385	0.0446546559342758	0.0334498027430146	1.54777780395155e-05	0.153439293417548	0.0144465160499012	-0.0742633817765345	0.0233028810049746	0.122615258145159	0.031585110811439	0.161220465932281	-0.0946599625516627	0.105815935459001	-0.160891824945655	0.0443439838308388];
% u_start = matUnstackH(u_start, 3);
u_start = stiefel_normalize(x, u_start);
[lambda_pim, v_pim] = pim_function(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigencheck_hessian(lambda_pim, v_pim, rhess_fun_han);

% lambda_pim = 960.3295;


if lambda_pim>0
    fprintf("lambda_pim %g\n", lambda_pim);
    
    mu = 1.1 * lambda_pim;

    rhess_shifted_fun_han = ...
        @(u) rsom_rhess_rot_stiefel(x,u,problem_struct_next) - mu .* u;
            
    %run shifted power iteration
    u_start_second_iter = stiefel_randTangentNormVector(x);
%     u_start_second_iter = [0.00077113264116923	0.281115776682463	-0.0743646397394423	-0.0886305960228544	0.117585914274081	-0.024036719437404	-0.260326468379064	-0.241679536479156	0.367204202483716	0.115337781039779	-0.674335601512407	0.164965856163922	-0.2088154371294	-0.0925590822378856	0.0689279154615746;
% 0.0696547686561622	-0.453027162805502	0.00402923055437596	-0.0466548942502283	-0.222891855621468	-0.157355948077402	-0.365962149166859	0.0222392492104101	-0.263012256818376	0.0738842664287403	-0.0397131251545104	-0.0316991930723567	-0.0830368274906565	0.305265399829179	-0.21406782865266;
% 0.279082504672713	-0.00567051550483593	-0.454409656879214	-0.195016878494919	0.0939708493354671	-0.136603144185852	0.212030090722661	0.00295869739043896	0.118278062847468	0.593606390806393	0.0155741214803132	-0.371875564251702	-0.319769654049333	-0.0149603955528836	0.0160555514771061;
% -0.318063536657995	0.562599862367956	-0.0570196471637393	0.673595895279574	0.12035604615016	0.606120881565591	0.494654711901422	-0.222430052693084	-0.426931204483174	0.0575233158356444	0.0424296566778161	-0.0265020704287926	-0.271903580565714	0.71437620951476	-0.331921083995758];
%     u_start_second_iter = matUnstackH(u_start_second_iter, 3);
    u_start_second_iter = stiefel_normalize(x, u_start_second_iter);
    [lambda_pim_after_shift, v_pim_after_shift] = pim_function( ...
        rhess_shifted_fun_han, u_start_second_iter, stiefel_normalize_han, thresh);
    
    disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
        'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
    eigencheck_hessian(lambda_pim_after_shift, v_pim_after_shift, ...
        rhess_shifted_fun_han);

    disp('Checking Eigenvalue shift:')
    disp(['difference between (lambda_pim_after_shift+mu)*v_pim_after_shift ' ...
        'and H(v_pim_after_shift) should be in the order of the tolerance:'])
    eigencheck_hessian(lambda_pim_after_shift + mu, v_pim_after_shift, ...
        rhess_fun_han);
end

%%%
disp(['Checking if highest_norm_eigenval = lambda_pim_after_shift + mu' ...
    ' is an eigenval for initial function:'])
disp(['difference between highest_norm_eigenval*v_pim and H(v_pim) ' ...
    'should be in the order of the tolerance:'])
highest_norm_eigenval = lambda_pim_after_shift + mu;
eigencheck_hessian(highest_norm_eigenval, v_pim, rhess_fun_han);
%%% scaling eigenvalue
% if ~eigencheck_hessian(highest_norm_eigenval, v_pim, rhess_fun_han)
%     % scale_factor
%     fac_1 = remove_quasi_zeros(highest_norm_eigenval*v_pim(:));
%     hess_hne = rhess_fun_han(v_pim);
%     fac_2 = remove_quasi_zeros(hess_hne(:));
%     fac2_1 = fac_2 ./ fac_1;
%     fac2_1_nums = fac2_1(~isnan(fac2_1));
%     fac2_1_finite = fac2_1_nums(isfinite(fac2_1_nums));
%     scale_factor = mode(fac2_1_finite);
% 
%     disp("Not even after scaling eigenval?")
%     eigencheck_hessian(scale_factor * highest_norm_eigenval, v_pim, rhess_fun_han);
%     highest_norm_eigenval = scale_factor * highest_norm_eigenval;
% end


%Preparing linesearch
nrs_next = problem_struct_next.sz(1);
d = problem_struct_next.sz(2);
N = problem_struct_next.sz(3);
step2.M = stiefelfactory(nrs_next,d,N);
step2.cost = @(x) rsom_cost_rot_stiefel(x, problem_struct_next);
step2.egrad = @(x) rsom_egrad_rot_stiefel(x, problem_struct_next);
step2.grad = @(x) rsom_rgrad_rot_stiefel(x, problem_struct_next);
step2.ehess = @(x,u) rsom_ehess_rot_stiefel(x,u, problem_struct_next);
step2.hess = @(x,u) rsom_rhess_rot_stiefel(x,u, problem_struct_next);


alphas = linspace(-0.01,0.01,501); %-0.2:0.01:0.2;
plot_vals = zeros(size(alphas));
plot_vals_taylor = zeros(size(alphas));
for ii = 1:length(alphas)
    x_retr_ii = step2.M.retr(x, v_pim_after_shift, alphas(ii));
%     disp("Is x_retr_ii on Stiefel? (Taylor)")
%     disp(check_is_on_stiefel(x_retr_ii));
%     disp([matStack(x), matStack(x_retr_ii)])
    plot_vals(ii) = step2.cost(x_retr_ii);
    %Note: gradient is zero
    plot_vals_taylor(ii) = step2.cost(x)+...
        alphas(ii)^2/2*sum(stiefel_metric(x,v_pim_after_shift,step2.hess(x,v_pim_after_shift),'canonical'));
end

plot(alphas, plot_vals,'b')
hold on
plot(alphas,plot_vals_taylor,'k.');
hold off


% alpha = min(lambdas_moved) + lambdas_max;
% alpha_linesearch = 10; %TODO: set this correctly
% SDPLRval = 10; %TODO: set this correctly 

disp("Now performing linesearch...");
%Note: first output param of linesearch() would be "stepsize"
[~, Y0] = linesearch_decrease(step2, ...
    x, -v_pim_after_shift, rsom_cost_rot_stiefel(x,problem_struct_next));

lambda_pim_out = lambda_pim_after_shift;
v_pim_out = v_pim_after_shift;


end %file function

