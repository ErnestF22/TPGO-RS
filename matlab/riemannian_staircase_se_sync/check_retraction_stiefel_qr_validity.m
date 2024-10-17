function b = check_retraction_stiefel_qr_validity(x,e)

% retraction_stiefel(x,e)
% x \in St(n,p)
% e \in T_x St(n,p)


p = size(x,1);
n = size(x,2);

%Check retraction validity:
%On Stiefel manifold, a retraction along direction $v \in T_{R0}St(n,p)$ 
% generates a curve $R(t)$ with $R(0)=R0$. 
%Per definition of retraction, the numerically computed derivative of 
% $R(t)$ for $t=0$ has to be equal to $v$ i.e., $\dot{R}(t)|_{t=0}=v$.

t = linspace(0.1, -0.1, 1001);

t_1 = floor(length(t));
t_2 = t_1+1;
h = t_2 - t_1;

% i_p = eye(p,p);
% snd_term = inv(sqrtm((i_p + e' * e)));
% rxe = (x + e)*snd_term;

rs = retraction_stiefel_qr(x,e);
rs_h = retraction_stiefel_qr(x + h, e);

num_der = (rs_h - rs) / h;

if (max(abs(e-num_der), [], "all")) < 1e-6
    b = boolean(1);
else
    b = boolean(0);
end


end % file function
