%% Assuming all vectors are defined in the natural metric coordinates, 
% IE:X,Y are vector fields on TSO(3)xSO(3) defined in the natural coordinates
clear all;

% Define metrics
M_nonnatural = randn(3);
M_nonnatural = M_nonnatural'*M_nonnatural;
[J,M_natural] = rotRef_SchurComplement(M_nonnatural);

% Create vector fields
x1Vec = randn(3,1); x2Vec = randn(3,1); x3Vec = randn(3,1);
y1Vec = randn(3,1); y2Vec = randn(3,1); y3Vec = randn(3,1);
X=@(R,U,RRef) [RRef*hat3(x3Vec);R*hat3(x1Vec);R*hat3(x2Vec)];
Y= @(R,U,RRef) [RRef*hat3(y3Vec);R*hat3(y1Vec);R*hat3(y2Vec)];

% Try parallel transporting vector at RRef to R, compute the nonnatural
% coordinates, then replace with original vector
x3AtR = @(R,RRef) rot_parallel(RRef,R,RRef*hat3(x3Vec),'torotation');
y3AtR = @(R,RRef) rot_parallel(RRef,R,RRef*hat3(y3Vec),'torotation');
X_R = @(R,U,RRef) [x3AtR(R,RRef);R*hat3(x1Vec);R*hat3(x2Vec)];
Y_R = @(R,U,RRef) [y3AtR(R,RRef);R*hat3(y1Vec);R*hat3(y2Vec)];

X_nonnatural_temp=@(R,U,RRef) kron(inv(J),eye(3))*X_R(R,U,RRef);
Y_nonnatural_temp=@(R,U,RRef) kron(inv(J),eye(3))*Y_R(R,U,RRef);
uVec = randn(3,1);
U=@(R) R*hat3(uVec);

% Define evaluation parameters
Reval = eye(3);%rot_randn;
RRefeval = rot_randn;
Ueval = U(Reval);

% Check that arclength is invariant
Z = [RRefeval;Reval;Ueval];
g_natural = rotRef_metric(Z,X(Reval,Ueval,RRefeval),...
    Y(Reval,Ueval,RRefeval),M_natural);

% Compute metric at the identity (natural coordinates)
g_natural_I = rotRef_metric([eye(3);eye(3);U(eye(3))],X(eye(3),U(eye(3)),eye(3)),...
    Y(eye(3),U(eye(3)),eye(3)),M_natural);

% Compute metric at the identity (nonnatural coordinates)
g_I_nonnatural = rotRef_metric_nonNatural([eye(3);eye(3);U(eye(3))],...
    kron(inv(J),eye(3))*X(eye(3),U(eye(3)),eye(3)),...
    kron(inv(J),eye(3))*Y(eye(3),U(eye(3)),eye(3)),M_nonnatural);
% replace the orignal vector for the SO(3) component
X_nonnatural = X_nonnatural_temp(Reval,Ueval,RRefeval);
X_nonnatural(1:3,:) = RRefeval*hat3(x3Vec);
Y_nonnatural = Y_nonnatural_temp(Reval,Ueval,RRefeval);
Y_nonnatural(1:3,:) = RRefeval*hat3(y3Vec);
g_nonnatural = rotRef_metric_nonNatural(Z,X_nonnatural,...
    Y_nonnatural,M_nonnatural);

fprintf('Assum Vector Fields are defined in natural coordinates,\n');
fprintf('natural, natural at Id, nonnatural, nonnatural at Id (should be equal):\n %0.3f, %0.3f, %0.3f, %0.3f\n',...
[g_natural g_natural_I g_nonnatural g_I_nonnatural]);

%% Assuming all vectors are defined in the nonnatural metric coordinates, 
% IE:X,Y are vector fields on TSO(3)xSO(3) defined in the nonnatural coordinates

% Compute using parallel transport
g_nonnatural_2 = rotRef_metric_nonNatural(Z,X(Reval,Ueval,RRefeval),...
    Y(Reval,Ueval,RRefeval), M_nonnatural);

% Use group property to move vector field to the Id
Z_I = [eye(3);eye(3);Reval'*Ueval];
X_I = rotRef_hat(eye(3),eye(3),rotRef_vee(Reval,RRefeval,X(Reval,Ueval,RRefeval)));
Y_I = rotRef_hat(eye(3),eye(3),rotRef_vee(Reval,RRefeval,Y(Reval,Ueval,RRefeval)));

% Compute at Id using left translation
g_I_nonnatural_2 = rotRef_metric_nonNatural(Z_I,X_I,Y_I,M_nonnatural);

% Compute at Id (natural method)
X_I_natural = kron(J,eye(3))*X_I;
Y_I_natural = kron(J,eye(3))*Y_I;
g_natural_I_2 = rotRef_metric(Z_I,X_I_natural,Y_I_natural,M_natural);

% Compute at R (natural method using parallel transport) (makes sense why
% this equals g_nonnatural_2 since they're both parallel transporting)
X_natural=kron(J,eye(3))*X_R(Reval,Ueval,RRefeval);
X_natural(1:3,:) = RRefeval*hat3(x3Vec);
Y_natural=kron(J,eye(3))*Y_R(Reval,Ueval,RRefeval);
Y_natural(1:3,:) = RRefeval*hat3(y3Vec);
g_natural_2 = rotRef_metric(Z,X_natural,Y_natural,M_natural);

fprintf('Assum Vector Fields are defined in nonnatural coordinates,\n');
fprintf('natural, natural at Id, nonnatural, nonnatural at Id (should be equal):\n %0.3f, %0.3f, %0.3f, %0.3f\n',...
[g_natural_2 g_natural_I_2 g_nonnatural_2 g_I_nonnatural_2]);
