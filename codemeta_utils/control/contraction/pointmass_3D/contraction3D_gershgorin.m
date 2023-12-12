function [ m, flag ] = contraction3D_gershgorin( eMin, eMax, kvIn, kdIn, betaIn )
%Gershgorin disc method to bound G [2 x 2] matrix
% INPUT:
% eMin - min eigenvalue of the Hessian of the cost function
% eMax - max eigenvalue of the Hessian of the cost function
% kvIn - scalar velocity gain
% kdIn - scalar position gain
% betaIn - scalar for contraction proof
% OUTPUT:
% m - [1 m12 m22] values
% flag - boolean where true means the system is feasible
%% Default Parameters
MAX_ITER = 100;
alpha = 10;
eps = 1e-3;
ERR_TOL = 1e-5;
flag = false; %set to false by default
% ksMax = 1;%eMax * kdIn; %using the max eigenvalue
ksMax = eMax*kdIn;
% ksMin = 0.1;%eMin * kdIn; %using the min eigenvalue
ksMin = eMin*kdIn;
ks=ksMin;
%% Optimize using gershgorin bounds
cvx_begin quiet
    variable m(2)
    Gmax = [2*m(1)+betaIn*m(2)-2*kvIn*m(2), betaIn*m(1)-ksMax*m(2)-kvIn*m(1)+1;...
        betaIn*m(1)-ksMax*m(2)-kvIn*m(1)+1, betaIn-2*ksMax*m(1)];
    Gmin = [2*m(1)+betaIn*m(2)-2*kvIn*m(2), betaIn*m(1)-ksMin*m(2)-kvIn*m(1)+1;...
        betaIn*m(1)-ksMin*m(2)-kvIn*m(1)+1, betaIn-2*ksMin*m(1)];
    minimize( max( [Gmax(1,1)-Gmax(1,2); Gmax(1,1)+Gmax(1,2); Gmax(2,2)-Gmax(2,1);Gmax(2,2)+Gmax(2,1);...
        Gmin(1,1)-Gmin(1,2); Gmin(1,1)+Gmin(1,2); Gmin(2,2)-Gmin(2,1);Gmin(2,2)+Gmin(2,1)] ) ) 
    subject to
        Gmax(1,1) <= 0%-abs(G(1,2)) %- 1e-3
        Gmax(2,2) <= 0%-abs(G(2,1)) %- 1e-3
        Gmin(1,1) <= 0%-abs(G(1,2)) %- 1e-3
        Gmin(2,2) <= 0%-abs(G(2,1)) %- 1e-3
        
        m(2)-m(1)^2 >= 0
cvx_end

% fprintf('cvx_optval: %f\n', cvx_optval);
if (~any(isnan(m)))
    eGmax = eig(Gmax);
% %     fprintf('eig(Gmax): %f, %f\n', eGmax(1), eGmax(2));
    eGmin = eig(Gmin);
%     fprintf('eig(Gmin): %f, %f\n', eGmax(1), eGmin(2));
    M = [1 m(1);m(1) m(2)];
    eM = eig(M);
%     fprintf('eig(M): %f, %f\n', eM(1), eM(2));
    if ( all(eGmax <= 0) && all(eig(M) > 0) && all(eGmin <= 0))
%         fprintf('System is feasible\n');
        flag = true;       
        %Check if found m's work for all eigenvalues in [eMin, eMax]
        %Since we've bounded G(Smin), G(Smax) to be less than 0, all S in
        %[Smin, Smax] will always be less than 0
%         contraction3D_checkGershgorin(eMin, eMax, kvIn, kdIn, betaIn, m); 
    end
end

end

