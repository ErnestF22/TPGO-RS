% Check if we can bound gershgorin disc radii by max eigenvalue of the
% non-diagonalizing matrix. IE if M = [A B;B D] where A=A', D=D' then
% S = [V 0;0 W] where V\A*V = diag(eig(A)) and W\D*W = diag(eig(D)). We
% want to see if gers. disc A+abs(B) <= max(eig(A)) + max((eig(B)).
% Proof:
% S\M*S = [V\A*V V\B*W;W\B*V W\D*W]
% Disc (A,B) = max(eig(A)) + abs(V\B*W)
%   abs(V\B*W)_{row 1} = ||V_{row 1}\B*W||_{1 norm}
%               ...
%               ...
%                      <= max(eig(B))
%  => max(eig(A)) + abs(V\B*W) <= max(eig(A)) + max(eig(B))

close all; clear all; clc;

for i = 1:1e6
    % Define matrices
    A = randn(3,3); A=A'*A; % A is symmetric
    D = randn(3,3); D=D'*D; % D is symmetric
    B = randn(3,3); B=B'*B; % B is arbit. symm matric

    % Compute diagonalizing matrices
    [V,A_eig] = eig(A);
    [W,D_eig] = eig(D);
    S = [V zeros(3);zeros(3) W];

    % Compute all gersh. disc centered at A
    newOffDiagMat = V\B*V;
    B_row1_abs = norm(newOffDiagMat(1,:),1);
    B_row2_abs = norm(newOffDiagMat(2,:),1);
    B_row3_abs = norm(newOffDiagMat(3,:),1);

    % Compute using bounds
    [B_vec,B_eig] = eig(B);
    B_radii_bound = sqrt(3)*max(abs(B_eig(:)));

    % Compare results
    fprintf('Analytical Eig,\t Bound,\t diff(should be positive)\n')
    results=[B_row1_abs, B_radii_bound, B_radii_bound-B_row1_abs;...
        B_row2_abs, B_radii_bound, B_radii_bound-B_row2_abs;...
        B_row3_abs, B_radii_bound, B_radii_bound-B_row3_abs]
    
    if any(results(:,3)<0)
        error('Bound does not work')
    end
end




        