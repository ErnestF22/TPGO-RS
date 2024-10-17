function POCdegree2Align
% From work with Ernesto
% See equations (41)-(45) in our working document
flagRandomizeQx=true;
flagPerturbMeasurements=false;

resetRands()

% Generate relative translations (with the last row already set to zero)
Tij=randn(3,2);
Ri=randrot_som(3);
Ritilde=[Ri;zeros(1,3)];
Tijtilde=Ritilde*Tij;
if flagPerturbMeasurements
    sigma=1e-3;
    Tijtilde=Tijtilde+sigma*randn(size(Tijtilde))
end

% Test computation of Qx
% It should align the last rows of Tijtilde to zero
Qx=align2d_EoF(Tijtilde);
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
disp([Ritilde*Tij Qx'*Qb*Qx*Tijtilde Tijtilde])
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
RbEst=procrustesRb_EoF(c,QLastRigh');
disp('Estimated and actual Rb should match up to a sign')
disp([RbEst Rb])

RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;

disp('One of the two residuals should be equal to zero')
disp(norm(RitildeEst1-Ritilde,'fro'))
disp(norm(RitildeEst2-Ritilde,'fro'))

% Repeat the check with the function that summarizes all the steps above
[RitildeEst1,RitildeEst2]=recoverRitilde_EoF(Ritilde2,Tijtilde);
disp('Again, one of the two residuals should be equal to zero')
disp(norm(RitildeEst1-Ritilde,'fro'))
disp(norm(RitildeEst2-Ritilde,'fro'))

function Qx=align2d_EoF(v)
Q=fliplr(orthComplement(v));
Qx=flipud(orthCompleteBasis(Q)');

function RbEst=procrustesRb_EoF(c,q)
[U,~,V]=svd(c*q');
RbEst=U*diag([1 det(U*V')])*V';

function [RitildeEst1,RitildeEst2]=recoverRitilde_EoF(Ritilde2,Tijtilde)
Qx=align2d_EoF(Tijtilde);
QxRitilde2Bot=Qx(3:4,:)*Ritilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRigh=Qx(3:4,4)';

RbEst=procrustesRb_EoF(c,QLastRigh');
RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
