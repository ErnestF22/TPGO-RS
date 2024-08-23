function test_POCdegree2Align
% From work with Ernesto
% See equations (41)-(45) in our working document
flagRandomizeQx=true;
flagPerturbMeasurements=false;

resetRands()

% Generate relative translations (with the last row already set to zero)
N = 6;
load('Qbnn_data/testdata.mat','edges')
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')


load('Qbnn_data/R_gt.mat','R_gt')
% load('Qbnn_data/testdata.mat','T')
load('Qbnn_data/testdata.mat','Tijs')

% nrs = size(R, 1);
% p = nrs;

R = R_gt;

d = size(R, 2);

% setup node degrees
low_deg = 2;
% node_degrees = [2,2,4,3,3,4];
% nodes_high_deg = node_degrees == 3 | node_degrees == 4;
% nodes_low_deg = ~nodes_high_deg;

%% Running procedure on node ii = 1

ii = 1;


%computing Tij12
Tij1j2 = zeros(d,low_deg);
num_js_found = 0;
for e = 1:size(edges)
    e_i = edges(e,1);
    %         e_j = edges(e,2);
    if (e_i == ii)
        num_js_found = num_js_found + 1;
        Tij1j2(:,num_js_found) = Tijs(:,e);
    end
end

% Tij=randn(3,2);
% Ri=randrot_som(3);
% T_edges = make_T_edges(T, edges);
% Ritilde=[Ri;zeros(1,3)];
Ritilde = [R(:,:,ii);zeros(1,3)];
Tijtilde=Ritilde*Tij1j2;
if flagPerturbMeasurements
    sigma=1e-3;
    Tijtilde=Tijtilde+sigma*randn(size(Tijtilde))
end

% Test computation of Qx
% It should align the last rows of Tijtilde to zero
Qx=align2d(Tijtilde);
if flagRandomizeQx
    % randomly combine the last two rows of Qx
    Rx=randrot_som(2);
    Qx=blkdiag(eye(2),Rx)*Qx;
end
disp(Qx*Tijtilde)

% Lemma: the top right 2x1 block of Qx is zero
% Proof: The last two rows of Qx need to span the orthogonal complement of
% span(Tijtilde). The top two rows of Qx need to span the ortogonal
% complement of the last two rows, i.e., they need to span the ortogonal
% complement of the orthogonal complement of span(Tijtilde), i.e.,
% span(Tijtilde) itself. However, by construction we have recovered
% Tijtilde to have zeros in the last row. Hence, also the top two rows of
% Qx need to have zeros in the last part (otherwise they would span a space
% different than span(Tijtilde)
% QED

% Generate a random Rb and the associated Qb
Rb=randrot_som(2);
Qb=blkdiag(eye(2),Rb);

% Check that Qx'*Qb*Qx leaves Tijtilde invariant
% And that the new Ritilde generated from the ambiguity
% still gives the same measurements.
% Ritilde2 is a 
Ritilde2=Qx'*Qb*Qx*Ritilde;
disp('The three sets of two columns each should be the same')
disp([Ritilde*Tij1j2 Qx'*Qb*Qx*Tijtilde Tijtilde])
% but Ritilde2 does not have the last row to zero
disp(Ritilde2)

% We want to find Rb, and Rb has the following properties
PLast=[0 0 0 1];
QLast=PLast*Qx';
disp(QLast*blkdiag(eye(2),Rb')*Qx*Ritilde2)

QxRitilde2Top=Qx(1:2,:)*Ritilde2;
QxRitilde2Bot=Qx(3:4,:)*Ritilde2;
QLastLeft=PLast*Qx(1:2,:)';
QLastRigh=PLast*Qx(3:4,:)';

% Lemma: QLastLeft is zero
% Proof: After taking the transpose, the 2x1 zero block in Qx mentioned in the
% previous Lemma becomes a 1x2 zero block in the bottom right. PLast*Qx'
% contains such block in the left part. Hence the claim.
% QED
disp('QLastLeft Should be zero')
disp(QLastLeft)

% Expanding the block diagonal matrix, we want to find Rb' such that the
% following is zero
disp('The expressions below should be a zero 1x3 block')
disp(QLastRigh*Rb'*QxRitilde2Bot+QLastLeft*QxRitilde2Top)
% However, using the previous Lemma, this can be simplified to
disp(QLastRigh*Rb'*QxRitilde2Bot)
% Lemma: QxRitilde2Bot (which is 2x3) has rank 1
% Proof: From the previous expression, Rb*QLastRigh' provides the
% coefficient of the linear combination of the rows of QxRitilde2Bot that
% gives a zero vector, showing that those two are not linearly independent
% QED

% Recover the coefficients of the linear combination from the proof above
% using SVD (up to a +/- sign!)
[U,S,~]=svd(QxRitilde2Bot,'econ');
disp('Last sval should be zero')
disp(S(2,2))

c=U(:,2);
disp('c and Rb*QLastRigh'' should match up to a +/- sign')
disp([c Rb*QLastRigh'])

% Estimate Rb by Procrustes alignment (without translation)
RbEst=procrustesRb(c,QLastRigh');
disp('Estimated and actual Rb should match up to a sign')
disp([RbEst Rb])

RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;

disp('One of the two residuals should be equal to zero')
disp(norm(RitildeEst1-Ritilde,'fro'))
disp(norm(RitildeEst2-Ritilde,'fro'))

% Repeat the check with the function that summarizes all the steps above
[RitildeEst1,RitildeEst2]=recoverRitilde(Ritilde2,Tijtilde);
disp('Again, one of the two residuals should be equal to zero')
disp(norm(RitildeEst1-Ritilde,'fro'))
disp(norm(RitildeEst2-Ritilde,'fro'))

disp([Qx' * blkdiag(eye(2),-RbEst') * Qx * Ritilde2, Ritilde]) % =\bmat{*_{3\times 3}\\0_{(p-3)\times 3}}
disp('')


function Qx=align2d(v)
Q=fliplr(orthComplement(v));
Qx=flipud(orthCompleteBasis(Q)');

function RbEst=procrustesRb(c,q)
[U,~,V]=svd(c*q');
RbEst=U*diag([1 det(U*V')])*V';

function [RitildeEst1,RitildeEst2]=recoverRitilde(Ritilde2,Tijtilde)
Qx=align2d(Tijtilde);
QxRitilde2Bot=Qx(3:4,:)*Ritilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRigh=Qx(3:4,4)';

RbEst=procrustesRb(c,QLastRigh');
RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
