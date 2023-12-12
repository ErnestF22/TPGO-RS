% Script to compare the results of the grid-bisection search algorithm in
% par_rotBundle_contractionTest.m. This script assumes the data has been
% loaded ( ie load('data_...mat') has been ran)

%% Find eigenvalue of Omega_D matrix
flag_printEigenValuesOfOmegaD = true;
syms w1Var w2Var w3Var m1Var m2Var m3Var kvVar real

wHat = [0 -w3Var w2Var;w3Var 0 -w1Var;-w2Var w1Var 0];

O11 = m2Var/4*wHat^2;
O21 = (m3Var*kvVar-m2Var)/4*wHat+m3Var/8*wHat^2;
O12 = O21';
O22 = zeros(3);

OmegaD = [O11 O12;O21 O22];
if flag_printEigenValuesOfOmegaD
    pretty(simplify(eig(OmegaD)));
end
    

% Replace OmegaD with real values
wMaxLinear = 40;
m1 = m_contract(1);m2 = m_contract(2); m3 = m_contract(4);
AA = vpa(subs(OmegaD,[m1Var;m2Var;m3Var;w1Var;w2Var;w3Var;kvVar],...
    [m1;m2;m3;wMaxLinear*ones(3,1);kv]));
E_OmegaD=eig(AA)

%% Find bounds P matrix
Lambda = diag([1,maxD/2*cot(maxD/2),maxD/2*cot(maxD/2)]);
PD11 = -m2*kd*Lambda + m1*beta*eye(3);
PD21 = -m3*kd/2*Lambda+(m2*beta+(m1-m2*kv)/2)*eye(3);
PD12 = PD21;
PD22 = (m3*beta+m2-m3*kv)*eye(3);
PD = [PD11 PD12;PD21 PD22];
E_PD=eig(PD)

% Find the actual contraction matrix
U = @(R) R*hat3(wMaxLinear*ones(3,1));
M = rotBundle_contractionMat(U,beta,kd,kv,[m1;m2;m3]);
maxE_contractionM = max(eig(M(R(:,:,1))+M(R(:,:,1))')/2)