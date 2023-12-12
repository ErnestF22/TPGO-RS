% Plot the results for the CDC 2019 paper "Geometric Attitude Control via Contraction on Manifold with
%   Automatic Gain Selection"
function [] = plotCDC2019(file2load)

if exist(file2load, 'file')
    data=load(file2load);
else
    data=load('data_contractionOpt_CDC2019.mat');
end

% Set figure formatting
setFigFontSize(8);
axisLimit = [-0.5 1];
stdFigSize = [3 1.25];

% Create figure for the grid-bisection search results
% figure
plot_surfPlot_contractionGains(data);
set(gca,'TickLabelInterpreter','latex')
savefigure('figGridResults_CDC2019','epsc',[3 3]);

% Plot the distance error on SO(3)
figure
normR = zeros(length(data.t),1);
for i = 1:length(data.t)
    normR(i) = norm(rot_vee(eye(3),rot3_log(data.R(:,:,i))));
end
plot(data.t,normR,'LineWidth',3);
% xlabel('Time (s)');
% ylabel('||(log_{I}R)^{v}||')
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{I_{3}}R)^{\vee}\|$','Interpreter','latex')
savefigure('figNormR_CDC2019','epsc',stdFigSize);

% Plot the norm(w) states
figure
normW = vecnorm(data.w);
plot(data.t,normW,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
% xlabel('Time (s)');
ylabel('$\|\omega\|$','Interpreter','latex')
savefigure('figNormV_CDC2019','epsc',stdFigSize);

% Plot max eigenvalues of contraction matrix
figure
m = [data.m_contract(1,1);data.m_contract(1,2);data.m_contract(2,2)];
symm = @(A) (A+A')/2; % return symmetric matric
for i = 1:length(data.t)
    U=@(R) R*hat3(data.w(:,i));
    M = rotBundle_contractionMat(U,data.beta,data.kd,data.kv,m);
    maxE(i) = max(eig(symm(M(data.R(:,:,i)))));
end 
plot(data.t,maxE,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
% xlabel('Time (s)');
ylabel('$\lambda_{max}(\mathcal{M})$','Interpreter','latex')
savefigure('figEig_CDC2019','epsc',stdFigSize);
end