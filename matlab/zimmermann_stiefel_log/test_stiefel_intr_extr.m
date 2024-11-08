clear;
clc;
close all;

x = make_rand_stiefel_3d_array(4,3,1);

[u,nu] = housegen(x);
disp('u');
disp(u);
disp('nu');
disp(nu);

[U,R] = hqrd(x);
disp('U');
disp(U);
disp('R');
disp(R);


%%
function [u,nu] = housegen(x)
    % [u,nu] = housegen(x)
    % Generate Householder reflection.
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    % [u,nu] = housegen(x).
    % H = I - uu' with Hx = -+ nu e_1
    %    returns nu = norm(x).
    u = x;
    nu = norm(x);
    if nu == 0
        u(1) = sqrt(2);
        return
    end
    u = x/nu;
    if u(1) >= 0
        u(1) = u(1) + 1;
        nu = -nu;
    else
        u(1) = u(1) - 1;
    end
    u = u/sqrt(abs(u(1)));
end
%%
function [U,R] = hqrd(X)  
    % Householder triangularization.  [U,R] = hqrd(X);
    % Generators of Householder reflections stored in U.
    % H_k = I - U(:,k)*U(:,k)'.
    % prod(H_m ... H_1)X = [R; 0]
    % where m = min(size(X))
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    [n,p] = size(X);
    U = zeros(size(X));
    m = min(n,p);
    R = zeros(m,m);
    for k = 1:min(n,p)
        [U(k:n,k),R(k,k)] = housegen(X(k:n,k));
        v = U(k:n,k)'*X(k:n,k+1:p);
        X(k:n,k+1:p) = X(k:n,k+1:p) - U(k:n,k)*v;
        R(k,k+1:p) = X(k,k+1:p);
    end
end

%%
function res_out = ObtainIntrHHR(x, etax, result)
    u = housegen(x);

    %HHRMTP
    % Return Q * Element or Q^T * Element or Element * Q or Element * Q^T, where T denote conjugate transpose
    % Q is H_1 H_2 .. H_k, and H_i is the householder reflector, such as that computed by HHRDecom.
    % lapack functions *unmqr_ or *ormqr_ are used.
    %  Trans: N, Side: L : Q * Element
    %  Trans: T, Side: L : Q^T * Element
    %  Trans: C, Side: L : Q^* * Element
    %  Trans: N, Side: R : Element * Q
    %  Trans: T, Side: R : Element * Q^T
    %  Trans: C, Side: R : Element * Q^*
         
    % tmp = HHRMtp()



end

% 
% Vector &Stiefel::ObtainIntrHHR(const Variable &x, const Vector &etax, Vector *result) const
% {
%     x.HHRDecom();
% 
%     tmp = HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::T, GLOBAL::L);
% 
%     for (integer i = 0; i < p; i++)
%     {
%         if (HHRptr[i + n * i] < 0)
%             scal_(&p, &GLOBAL::DNONE, tmpptr + i, &n);
%     }
% 
%     realdp *resultptr = result->ObtainWriteEntireData();
% 
%     realdp r2 = static_cast<realdp>(sqrt(2.0));
%     integer idx = 0;
%     for (integer i = 0; i < p; i++)
%     {
%         for (integer j = i + 1; j < p; j++)
%         {
%             resultptr[idx] = r2 * (tmpptr[j + i * n] - tmpptr[i + j * n]) / 2;
%             idx++;
%         }
%     }
% 
%     for (integer i = 0; i < p; i++)
%     {
%         for (integer j = p; j < n; j++)
%         {
%             resultptr[idx] = tmpptr[j + i * n];
%             idx++;
%         }
%     }
%     return *result;
% };
% 
% void Element::HHRDecom(void) const
%     {
%         if(FieldsExist("_HHR") && FieldsExist("_tau"))
%             return;
% 
%         if(!iscomplex)
%         {
%             integer m = size[0], n = size[1], minmn = (m < n) ? m : n;
%             Element HHR(*this), tau(minmn);
%             realdp *HHRptr = HHR.ObtainWritePartialData();
%             realdp *tauptr = tau.ObtainWriteEntireData();
%             integer *jpvt = new integer[n];
%             integer info;
%             integer lwork = -1;
%             realdp lworkopt;
%             for (integer i = 0; i < n; i++)
%                 jpvt[i] = i + 1;
%             /*  compute the space required in the geqp3 */
%             geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, &lworkopt, &lwork, &info);
%             lwork = static_cast<integer> (lworkopt);
%             realdp *work = new realdp[lwork];
%             /* QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
%             details: www.netlib.org/lapack/explore-html/db/de5/geqp3_8f.html */
%             geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, work, &lwork, &info);
%             if (info < 0)
%                 printf("Error in qr decomposition!\n");
% 
%             for (integer i = 0; i < n; i++)
%             {
%                 if (jpvt[i] != (i + 1))
%                     printf("Error in qr decomposition!\n");
%             }
%             AddToFields("_HHR", HHR);
%             AddToFields("_tau", tau);
%             delete[] jpvt;
%             delete[] work;
%             return;
%         }
% 
%         integer m = size[0] / 2, n = size[1], minmn = (m < n) ? m : n;
%         Element HHR(*this), tau(minmn, 1, 1, "complexRopt");
%         realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainWritePartialData();
%         realdpcomplex *tauptr = (realdpcomplex *) tau.ObtainWriteEntireData();
%         integer *jpvt = new integer[n];
%         integer info;
%         integer lwork = -1;
%         realdpcomplex lworkopt;
%         for (integer i = 0; i < n; i++)
%             jpvt[i] = i + 1;
%         realdp *rwork = new realdp[2 * n];
%         /* compute the space required in the geqp3 */
%         geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, &lworkopt, &lwork, rwork, &info);
%         lwork = static_cast<integer> (lworkopt.r);
%         realdpcomplex *work = new realdpcomplex[lwork];
%         /* QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
%         details: www.netlib.org/lapack/explore-html/db/de5/geqp3_8f.html */
%         geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, work, &lwork, rwork, &info);
%         if (info < 0)
%             printf("Error in qr decomposition!\n");
% 
%         for (integer i = 0; i < n; i++)
%         {
%             if (jpvt[i] != (i + 1))
%                 printf("Error in qr decomposition!\n");
%         }
%         AddToFields("_HHR", HHR);
%         AddToFields("_tau", tau);
%         delete[] jpvt;
%         delete[] work;
%         delete[] rwork;
%     };
% 
% 
%     Element Element::HHRMtp(Element HHR, Element tau, char *Trans, char *Side) const
%     {
%         assert(iscomplex == HHR.Getiscomplex() && iscomplex == tau.Getiscomplex());
% 
%         if(!iscomplex)
%         {
%             integer m = size[0], n = size[1], r = HHR.Getsize()[0], k = HHR.Getsize()[1], minrk = (r < k) ? r : k;
%             Element result(*this);
%             realdp *resultptr = result.ObtainWritePartialData();
%             const realdp *HHRptr = HHR.ObtainReadData();
%             const realdp *tauptr = tau.ObtainReadData();
% 
%             integer info;
%             integer lwork = -1;
%             realdp lworkopt;
%             /* compute the size of space required in the ormqr */
%             ormqr_(Side, Trans, &m, &n, &minrk, const_cast<realdp *> (HHRptr), &r, const_cast<realdp *> (tauptr), resultptr, &m, &lworkopt, &lwork, &info);
%             lwork = static_cast<integer> (lworkopt);
%             realdp *work = new realdp[lwork];
% 
%             ormqr_(Side, Trans, &m, &n, &minrk, const_cast<realdp *> (HHRptr), &r, const_cast<realdp *> (tauptr), resultptr, &m, work, &lwork, &info);
%             delete [] work;
% 
%             return result;
%         }
% 
%         integer m = size[0] / 2, n = size[1], r = HHR.Getsize()[0] / 2, k = HHR.Getsize()[1], minrk = (r < k) ? r : k;
%         Element result(*this);
%         realdpcomplex *resultptr = (realdpcomplex *) result.ObtainWritePartialData();
%         const realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainReadData();
%         const realdpcomplex *tauptr = (realdpcomplex *) tau.ObtainReadData();
% 
%         integer info;
%         integer lwork = -1;
%         realdpcomplex lworkopt;
% 
%         unmqr_(Side, Trans, &m, &n, &minrk, const_cast<realdpcomplex *> (HHRptr), &r, const_cast<realdpcomplex *> (tauptr), resultptr, &m, &lworkopt, &lwork, &info);
%         lwork = static_cast<integer> (lworkopt.r);
%         realdpcomplex *work = new realdpcomplex[lwork];
% 
%         unmqr_(Side, Trans, &m, &n, &minrk, const_cast<realdpcomplex *> (HHRptr), &r, const_cast<realdpcomplex *> (tauptr), resultptr, &m, work, &lwork, &info);
%         delete [] work;
% 
%         return result;
%     };
