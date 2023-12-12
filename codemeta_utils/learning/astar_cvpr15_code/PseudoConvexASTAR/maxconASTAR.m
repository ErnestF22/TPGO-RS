
% ##########################################################################
%% G L O B A L L Y   O P T I M A L   C O N S E N S U S   M A X I M I S A T I O N
%% This package contains the source code which implements optimal Consensus 
% Maximisation proposed in
% T.J. Chin, P. Purkait, A. Eriksson and D. Suter
% Efficient Globally Optimal Consensus Maximisation with Tree Search, 
% In Proceedings of the IEEE Conference on Computer Vision and Pattern 
% Recognition (CVPR), June 2015, Boston
% 
% Copyright (c) 2015 Pulak Purkait (pulak.purkait@adelaide.edu.au.)
% School of Computer Science, The University of Adelaide, Australia
% The Australian Center for Visual Technologies
% http://www.cs.adelaide.edu.au/directory/pulak.purkait
%% Please acknowledge the authors by citing the above paper in any academic 
%  publications that have made use of this package or part of it.
% ##########################################################################

function [pk, wk, vk,nnum,xnum] = maxconASTAR(A, b, c, d, th)
% BFS with Matousek for maximum consensus.

% fmincon options
optmincon = optimoptions(@fmincon, 'MaxIter', 500, 'Algorithm','active-set'); 
optmincon = optimoptions(optmincon,'Display',  'off', 'Diagnostics', 'off');%'iter'); 
% optmincon = optimoptions(optmincon,'GradObj','on','GradConstr','on');%, 'Hessian', 'user-supplied', 'HessFcn', @hessinterior);%, 'TolX', 0, 'TolFun', 0);

noimages = numel(d); 
p = size(A, 2); 
x0 = [randn(p, 1); 0]; 
I0 = [1:noimages]; 

xnum = 1; 
% x0 = rand(2, 1); 
[node.p, node.w, node.b] = solve_fmincon_minimax(A, b, c, d, optmincon, x0); 
x0 = [node.p; node.w]; 
 
dictHstc = java.util.Hashtable; 
node.v = [];   % Violation set.

[node.h, w, p, ~, ~, xnum,dictHstc] = heuristic(A, b, c, d, optmincon, I0, [], [], x0, th, xnum,dictHstc);     % Heuristic value.(X,y, H, v, x0, cyc, th)
[U_cnt, v] = compute_upper(A, b, c, d, p, th); 

ps = numel(p); 

node.f = length(node.v)+node.h; % Evaluation value.
nnum = 0;
if (node.w<=th) 
    pk = node.p;
    wk = node.w;
    vk = node.v; 
    return;
end
if  node.h >= U_cnt
    pk = p;
    wk = w;
    vk = v; 
    return;
end


node.solved = 0; 
minOUT = U_cnt; 
q = node; 

priof = node.f;
priow = node.w; 
% wprio = node.w; 
% List of checked bases.
dictS = java.util.Hashtable; 
% dict = zeros(1, 10^3*N); 

while 1

    [mn,~] = min(priof, [], 2);
    id = find(priof == mn); 
    [~, mp] = min(priow(id));
    mp = id(mp); 
 
    node = q(mp);
    if(node.solved) || node.w<=th  % Reached a goal state.
        pk = node.p;
        wk = node.w;
        vk = node.v; 
        return; 
    end
    q(mp) = [];
    
    disp([numel(node.v), mn, minOUT]); 
%     disp([mp, priof]); 
    priof(mp) = []; 
    priow(mp) = []; 
    nnum = nnum+1;
            
    % Generate child nodes.
    for i=1:length(node.b)
        xnum = xnum+1;

%         mat2str(node.v); 
%         num2str(node.v); 
%         child.v = [node.v node.b(i)];
        child.v = sort([node.v node.b(i)]);
        entry = dictS.get(sprintf('%4ld', child.v)); 
        if numel(entry)
            continue; 
        end
        
        H = 1:noimages;
        H(child.v) = []; 

        
        entry = dictHstc.get(sprintf('%4ld', child.v));
        if numel(entry)
            child.h = entry(ps+2); 
            w = entry(ps+1); 
            p = entry(1:ps); 

        else
            Ih = [2*(H-1)+1; 2*H]; 
            [child.p, child.w, bs] = solve_fmincon_minimax(A(Ih(:), :), b(:, H), c(:, H), d(H), optmincon, x0);
            child.b = H(bs);
            [child.h, w, p, ~, ~, xnum,dictHstc] = heuristic(A, b, c, d, optmincon, H, [], child.v, x0, th, xnum,dictHstc);     % Heuristic value.(X,y, H, v, x0, cyc, th)
            
        end
        
        child.h = max(child.h, node.h-1); 
        child.f = numel(child.v)+child.h;

        dictS.put(sprintf('%4ld', child.v), [1, 1]); 
        [UN_cnt, v] = compute_upper(A(Ih(:), :), b(:, H), c(:, H), d(H), p, th); 
        if  minOUT > UN_cnt + numel(child.v)
            minOUT = UN_cnt + numel(child.v); 
        end

        if child.h >= UN_cnt && child.f <= minOUT %|| child.h < 2
            % Reached a to the solution; 
            child.solved = 1; 
            child.p = p; 
            child.w = w; 
            child.v = sort([child.v, H(v)]); 
        else
            child.solved = 0; 
            if child.w<=th
                % Reached a goal state.
                pk = node.p;
                wk = node.w;
                vk = child.v; 
                return;
            end
        end
            idx = priof > minOUT; 
            q(idx) = []; 
            priof(idx) = []; 
            priow(idx) = []; 
            
            q = [ q, child ]; %#ok<AGROW>
            priof = [priof, child.f]; 
            priow = [priow, child.w];  
        
    end
end

end


function [ cyc, w, x0_f, chsub_old, sub_f, xnum, dictHstc] = ...
    heuristic(A, b, c, d, optmincon, H, chsub_old, v, x0, th, xnum, dictHstc)

Ih = [2*(H-1)+1; 2*H]; 
[x0, val, chsub] = solve_fmincon_minimax(A(Ih(:), :), b(:, H), c(:, H), d(H), optmincon, x0);
xnum = xnum + 1; 

% r = X(chsub_old, :)*x0 - y(chsub_old, :); 
% M = find(abs(r) < val);
% H = [H, chsub_old(M)]; 
% chsub_old(M) = []; 

if val>th
    % Feasible solution not found.
    H1 = H; 
    H1(chsub) = [];
    
     entry = dictHstc.get(sprintf('%4ld', sort([v, H(chsub)])));
    if numel(entry)
        p = numel(x0); 
        cyc = entry(p+2); 
        w = entry(p+1); 
        x0_f = entry(1:p); 
        sub_f = entry(p+3:end)'; 
        chsub_new = chsub_old; 
        disp('hash'); 
    else  
        [ cyc, w, x0_f, chsub_new, sub_f, xnum, dictHstc] = ...
            heuristic(A, b, c, d, optmincon, H1, H(chsub), sort([v, H(chsub)]), [x0; val], th, xnum, dictHstc); 
    end
else
    cyc = 0;
    x0_f = x0; 
    w = val; 
    sub_f = H; 
    return; 
end
% disp(cyc); 
for i=1:numel(chsub_new)
    sub_f = [sub_f, chsub_new(i)]; 
    if numel(sub_f) <= numel(chsub)
          disp(numel(sub_f)); 
          cyc = cyc + 1;
          return; 
    end
    Ih = [2*(sub_f-1)+1; 2*sub_f]; 
    
    resn = compute_residuals(A(Ih(:), :), b(:, sub_f), c(:, sub_f), d(sub_f), x0);
    if (max(resn) - th > 0.0001) || min(resn) < 0
        [x1, valn, bs] = solve_fmincon_minimax(A(Ih(:), :), b(:, sub_f), c(:, sub_f), d(sub_f), optmincon, [x0; max(resn)]);
        xnum = xnum + 1; 
    else
        continue; 
    end
    
    if valn > th
        pts = sub_f(bs); 
        sub_f(bs) = []; 
        
        if numel(sub_f) <= numel(chsub)
            disp(numel(sub_f)); 
            cyc = cyc + 1;
            return; 
        end
        
        cyc = cyc + 1; 
        
    else
        x0_f = x1;
        w = valn; 
    end
        
end
if numel(v)
    dictHstc.put(sprintf('%4ld', v), [x0_f; w; cyc;sub_f']); 
end

end

function resn = compute_residuals(A, b, c, d, x)
    nbrimages = numel(d); 
    p = size(A, 2); 
    AxPb = reshape(A*x(1:p),2,nbrimages)+b; % Ax + b
    Sqrt_AxPb = sqrt(sum(AxPb.^2));                     % sqrt(Ax + b)
    CxPd = (x(1:p)'*c + d)+eps; 
    resn =  zeros(1, nbrimages); 
    id = abs(CxPd) > 0.0001; 
    resn(id) = Sqrt_AxPb(id)./CxPd(id);
    resn(~id) = max(resn);
end

function [cnt, M] = compute_upper(A, b, c, d, x, valn)
        nbrimages = numel(d); 
        p = size(A, 2); 
        AxPb = reshape(A*x(1:p),2,nbrimages)+b; % Ax + b
        Sqrt_AxPb = sqrt(sum(AxPb.^2));                     % sqrt(Ax + b)
        CxPd = (x(1:p)'*c + d)+eps; 
        resn =  zeros(1, nbrimages); 
        id = abs(CxPd) > 0.0001; 
        resn(id) = Sqrt_AxPb(id)./CxPd(id);
    %     resn(~id) = x(end); 
        resn(~id) = max(resn);
        M = find(resn > valn); 
        cnt = numel(M); 
end

