function [ output_args ] = contraction3D_checkGershgorin( eMin, eMax, kvIn, kdIn, betaIn, m )
% Check if the given m values works for all evalues in [eMin, eMax] in the
% G matrix
% INPUT:
% eMin - min eigenvalue of the Hessian of the cost function
% eMax - max eigenvalue of the Hessian of the cost function
% kvIn - scalar velocity gain
% kdIn - scalar position gain
% betaIn - scalar for contraction proof (max convergence rate)
% m - a [2x1] vector of [m12; m22] values

%% Default Parameters
stepSize = 1e-3;

%% Visual Plot of f(s)
% f(s) = max( G(1,1)[s] + abs(G(1,2)[s]), G(2,2)[s] - abs(G(2,1)[s]) )
figure
title('f(s) = max(G(1,1)+abs(G(1,2)), G(2,2)+abs(G(2,1)))');
hold on
for i = eMin:stepSize:eMax
    ks = i*kdIn;
    G = [2*m(1)+betaIn*m(2)-2*kvIn*m(2), betaIn*m(1)-ks*m(2)-kvIn*m(1)+1;...
        betaIn*m(1)-ks*m(2)-kvIn*m(1)+1, betaIn-2*ks*m(1)];
    f = max( [G(1,1) + abs(G(1,2)), G(2,2) + abs(G(2,1))] );
    if (f<0)
        plot(i, f, 'gx');
    else
        plot(i, f, 'rx');
    end
end

%% find min/max points
%There can be a maximum of 4 pivot points located at ks = [eMin*kd, eMax*kd,
%alpha, gamma];

%alpha is the ks value for which G(1,1) = G(2,2), point at which the lesser
%eigenvalue becomes the greater
alpha = (betaIn-2*m(1)-betaIn*m(2)+2*kvIn*m(2))/(2*m(1));
%gamma is the ks value for which |B*m12-ks*m22-kv*m12+1| = 0, switching
%point for the off diagonal term (only one exist since G is symmetric)
gamma = (betaIn*m(1)-kvIn*m(1)+1)/(m(2));

%create list of ks values to check for max of f
ks_list = [eMin*kdIn, eMax*kdIn];

%check if alpha is within the bounds of (eMin, eMax)
if (alpha < eMax && alpha > eMin)
    ks_list = [ks_list, alpha];
end

%check if gamma is within the bounds of (eMin, eMax)
if (gamma < eMax && gamma > eMin)
    ks_list = [ks_list, gamma];
end

%Calculate f for each pivot point
max_f_list = [];
for i = 1:length(ks_list)
    %Define G
    ks = ks_list(i);
    G = [2*m(1)+betaIn*m(2)-2*kvIn*m(2), betaIn*m(1)-ks*m(2)-kvIn*m(1)+1;...
        betaIn*m(1)-ks*m(2)-kvIn*m(1)+1, betaIn-2*ks*m(1)];
    %Check conditions to determine which eigenvalue is greatest
    if (m(1) > 0)
        if ( ks_list(i) > alpha )
            %G(1,1) + abs(G(1,2)) > G(2,2)+abs(G(2,1))
            max_f_list(i) = G(1,1)+abs(G(1,2));
        elseif ( ks_list(i) < alpha)
            %G(1,1) + abs(G(1,2)) < G(2,2)+abs(G(2,1))
            max_f_list(i) = G(2,2)+abs(G(2,1));
        else
            %Both eigenvalues are equal, use G(1,1)+abs(G(1,2))
            max_f_list(i) = G(1,1)+abs(G(1,2));
        end
    elseif (m(1) < 0)
        if ( ks_list(i) > alpha )
            %G(1,1) + abs(G(1,2)) < G(2,2)+abs(G(2,1))
            max_f_list(i) = G(2,2)+abs(G(2,1));
        elseif ( ks_list(i) < alpha)           
            %G(1,1) + abs(G(1,2)) > G(2,2)+abs(G(2,1))
            max_f_list(i) = G(1,1)+abs(G(1,2));
        else
            %Both eigenvalues are equal, use G(1,1)+abs(G(1,2))
            max_f_list(i) = G(1,1)+abs(G(1,2));
        end
    else
        error('m(1) = 0: not implemented');
    end
end

[maxF, idxF] = max(max_f_list);
fprintf('Max f = %f at ks = %f (eVal = %f)\n', maxF, ks_list(idxF), ks_list(idxF)/kdIn);
if (maxF < 0)
    fprintf('Max f < 0: True\n');
else
    fprintf('Max f < 0: False\n');
end
end

