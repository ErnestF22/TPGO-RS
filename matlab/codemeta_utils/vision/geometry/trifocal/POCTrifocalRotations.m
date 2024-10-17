function POCTrifocalRotations
data=trifocal_dataset();

nTask=7;

switch nTask
    case 1
        %check trifocal line constraint
        R21=data.R(:,:,1)*data.R(:,:,2)';
        R23=data.R(:,:,2)*data.R(:,:,3)';

        checkConstraint(R21,R23,data)
    case 2
        %direct quaternion polynomial embedding
        l=sym('l',[3 3]);
        q=sym('q',[4,2]);

        R12q=quat_toR(q(:,1));
        R23q=quat_toR(q(:,2));
        l1=l(:,1);
        l2=l(:,2);
        l3=l(:,3);
        p=collect(l1.'*R12q*hat(l2)*R23q*l3,q);
        [c,t]=coeffs(p,q);
        %matlabFunction(c,'file','trifocal_quatPolynomial_coeffAuto');
        %matlabFunction(t,'file','trifocal_quatPolynomial_embeddingAuto');
    case 3
        %double quadratic embedding
        l=sym('l',[3 3]);
        R=sym('R',[6 3]);
        rho=sym('rho');
        x=[R(:); rho];

        R21=R(1:3,:);
        R23=R(4:6,:);
        l1=l(:,1);
        l2=l(:,2);
        l3=l(:,3);
        disp('# Coefficients and embedding for data equation + rotation constraints')
        p1=collect(l1.'*R21*hat(l2)*R23*l3,R);
        p2a=R21.'*R21-rho*eye(3);
        p2b=R21*R21.'-rho*eye(3);
        p3a=R23.'*R23-rho*eye(3);
        p3b=R23*R23.'-rho*eye(3);
        p=[p1;p2a(:);p2b(:);p3a(:);p3b(:)];
        [c,t]=coeffsArray(p,x);
        Nt=length(t);
        symMatlabFunctionFile(c{1},'trifocal_line_dataConstraints','l');
        symMatlabFunctionFile(cat(1,c{2:end}),'trifocal_line_identityConstraints','');
        symMatlabFunctionFile(subs(t,rho,1)','trifocal_line_embedding','R');
        disp('# Converting kronecker square of embedding to strings')
        Ntt=Nt*(Nt+1)/2;
        tt=sym(zeros(Ntt,1));
        ttIdx=zeros(Ntt,2);
        ttStr=repmat(' ',Ntt,1);
        cnt=1;
        w=getTextWaitBar(Ntt);
        w(0)
        for it=1:Nt
            for jt=it:Nt
                tt(cnt)=t(it)*t(jt);
                ttIdx(cnt,:)=[it jt];
                s=char(tt(cnt));
                ttStr(cnt,1:length(s))=s;
                w(cnt)
                cnt=cnt+1;
            end
        end
        ttStr(ttStr==0)=' ';
%         idxMat=repmat(1:Nt,Nt,1);
%         ttIdx=[reshape(idxMat,[],1) reshape(idxMat',[],1)];
%         tt=kron(t,t);
%         Ntt=length(tt);
%         ttStr=sym2char(tt,'progressBar');
        disp('# Finding equality constraints')
        [ttStrSorted,idxStrTt]=sortrows(ttStr);
        idxEqTt=[];
        w=getTextWaitBar(Ntt);
        w(0)
        for itt=2:Ntt
            if strcmp(ttStrSorted(itt,:),ttStrSorted(itt-1,:))
                idxEqTt=[idxEqTt; idxStrTt([itt-1 itt])'];
            end
            w(itt)
        end
        NEquality=size(idxEqTt,1);
        disp([num2str(NEquality) ' equality constraints'])
        cEquality=repmat([1 -1],[1,1,NEquality]);
        idxEquality=zeros(2,2,NEquality);
        for iEquality=1:NEquality
            idxEquality(:,:,iEquality)=[ttIdx(idxEqTt(iEquality,1),:)' ttIdx(idxEqTt(iEquality,2),:)'];
        end
        %the coefficient for lambda_ab for the i-th equation is obtained as 
        %\sum_j c(1,j,i)*v_a(idx(1,j,i))*v_b(idx(2,j,i))
        disp('# Finding rotation constraints')
        Zpoly=[];
        for iIdentity=1:8
            switch iIdentity
                case 1
                    Ipolyi=reshape(expand(R21*R21.'*R23*R23.'),[],1);
                case 2
                    Ipolyi=reshape(expand(R21.'*R21*R23*R23.'),[],1);
                case 3
                    Ipolyi=reshape(expand(R21*R21.'*R23.'*R23),[],1);
                case 4
                    Ipolyi=reshape(expand(R21.'*R21*R23.'*R23),[],1);
                case 8
                    Ipolyi=reshape(expand(R21.'*R23*R23.'*R21),[],1);
%                 case 5
%                     Ipolyi=reshape(expand(R12*R23.'*R23*R12.'),[],1);
%                 case 6
%                     Ipolyi=reshape(expand(R12.'*R23*R23*R12.'),[],1);
%                 case 7
%                     Ipolyi=reshape(expand(R12*R23.'*R23.'*R12),[],1);
            end
            Zpoly=[Zpoly;Ipolyi([2 3 4 6 7 8])];
        end
        
        %assumes all constraints use only the available monomials
        NIdentity=length(Zpoly);
        cIdentity=ones(1,27,NIdentity);
        idxIdentity=zeros(2,27,NIdentity);
        for iIdentity=1:NIdentity
            [cZ,tZ]=coeffs(Zpoly(iIdentity),x);
            tZStr=sym2char(tZ);
            [~,idxIntersectionZ,idxIntersectionAll]=intersect(tZStr,ttStr,'rows');
            if length(idxIntersectionZ)<length(tZ);
                error('Identity constraint use monomials not present before')
            end
            idxIdentity(:,:,iIdentity)=ttIdx(idxIntersectionAll,:)';
        end
        
        save('trifocal_line_coefficients','cEquality','idxEquality','cIdentity','idxIdentity')
%         [cZ,tZ]=coeffsArray(Zpoly,x,'progressBar');
%         NIdentityEq=length(cZ);
%         NMonomialsEq=length(tZ);
%         tZStr=sym2char(tZ,'progressBar');
%         flagMonomialPresent=ismember(tZStr,ttStr,'rows');
%         if all(flagMonomialPresent)
%             disp('All monomials present')
%         else
%             error('Identity constraints use monomials not present before')
%         end
%         cZDoubleCell=cellfun(@double,cZ, 'UniformOutput', false);
%         cZDouble=cat(1,cZDoubleCell{:});
%         [strIntersection,idxIntersectionZ,idxIntersectionAll]=intersect(tZStr,ttStr,'rows');
% %         flagIdentityEquValid=all(~and(cZDouble==1,repmat(flagMonomialPresent',NIdentityEq,1)==0),2);
    case 4
        %check constraint functions
        R21=data.R(:,:,1)*data.R(:,:,2)';
        R23=data.R(:,:,2)*data.R(:,:,3)';
        t=trifocal_line_embedding([R21;R23]);
        NLines=size(data.l2d,2);
        M=zeros(NLines,length(t));
        for iLine=1:NLines
            M(iLine,:)=trifocal_line_dataConstraints(squeeze(data.l2d(:,iLine,:)));
        end
        M=[M; trifocal_line_identityConstraints];
        disp(max(abs(M*t)))
        [U,S,V]=svd(M);
        szM=size(M);
        kerDim=max(szM(2)-szM(1),0)+sum(diag(S)<1e-12);
        kerBasis=V(:,end-kerDim+1:end);
        lambda=kerBasis'*t;
        disp(max(abs(t-kerBasis*lambda)))
        
        tLambda=trifocal_line_lambda_embedding(lambda);
        load('trifocal_line_coefficients');
        AEquality=trifocal_line_lambda_buildEquations(cEquality,idxEquality,kerBasis);
        AIdentity=trifocal_line_lambda_buildEquations(cIdentity,idxIdentity,kerBasis);
        A=[AEquality;AIdentity];
        disp(max(abs(A*tLambda)))
    case 5
        %test to extract rotations with rectification of views using QR
        %decomposition
        l2d1=data.l2d(:,:,1);
        l2d2=data.l2d(:,:,2);
        l2d3=data.l2d(:,:,3);
        [Q1,l2d1]=flipqr(l2d1);
        [Q2,l2d2]=flipqr(l2d2);
        [Q3,l2d3]=flipqr(l2d3);
        l2d1=homogeneous(l2d1,3);
        l2d2=homogeneous(l2d2,3);
        l2d3=homogeneous(l2d3,3);
        R21=Q1'*data.R(:,:,1)*data.R(:,:,2)'*Q2;
        R23=Q2'*data.R(:,:,2)*data.R(:,:,3)'*Q3;

        NLines=size(data.l3d,3);
        for iLine=1:NLines
            l1=l2d1(:,iLine);
            l2=l2d2(:,iLine);
            l3=l2d3(:,iLine);
            disp(l1'*R21*hat(l2)*R23*l3)
        end
        
    case 6
        %constraints and their derivatives
        R21=data.R(:,:,2)*data.R(:,:,1)';
        R23=data.R(:,:,2)*data.R(:,:,3)';

        l1=data.l2d(:,1,1);
        l2=data.l2d(:,1,2);
        l3=data.l2d(:,1,3);

        disp(l1'*R21'*hat(l2)*R23*l3)
        disp(l3'*R23'*hat(l2)*R21*l1)
        g=@(R21,R23) [hat(l1)*R21'*hat(l2)*R23*l3;hat(l3)'*R23'*hat(l2)*R21*l1];
        g0=g(R21,R23);
        v0=cnormalize(g0);
        v=(eye(6)-v0*v0')*randn(6,1);
        v1=v(1:3);
        v2=v(4:6);
        R21t=rot_geodFun(R21,rot_hat(R21,v1));
        R23t=rot_geodFun(R23,rot_hat(R23,v2));
        
        f=@(t) l1'*R21t(t)'*hat(l2)*R23t(t)*l3;
        df=@(t) v'*g(R21t(t),R23t(t));
        
        check_der(f,df)
    case 7
        %all constraints as SDP relaxation
        l=sym('l',[3 3]);
        r=sym('R',[18 1]);
        rho=1;%sym('rho');
        %x=[r; rho];

        R21=reshape(r(1:9),[3 3]);
        R23=reshape(r(10:18), [3 3]);
        l1=l(:,1);
        l2=l(:,2);
        l3=l(:,3);
        disp('# Coefficients and embedding for data equation + rotation constraints')
        p1=collect(l1.'*R21*hat(l2)*R23*l3,r);
        p2a=R21.'*R21-rho*eye(3);
        p2b=R21*R21.'-rho*eye(3);
        p3a=R23.'*R23-rho*eye(3);
        p3b=R23*R23.'-rho*eye(3);
        p=[p1;p2a(:);p2b(:);p3a(:);p3b(:)];
        [c,t]=coeffsArray(p,r);
        %disp('# Writing files to generate coefficients')
        %symMatlabFunctionFile(c{1},'trifocal_line_preSDP_dataConstraints','l');
        %symMatlabFunctionFile(cat(1,c{2:end}),'trifocal_line_preSDP_identityConstraints','');
        disp('# Computing maps from vector to matrix for the coefficients')
        Nt=length(t);
        M=zeros(Nt,18*18);
        m=0;
        fid=fopen('trifocal_line_preSDP_mappingMatrix.m','wt');
        fprintf(fid,'function [M,m]=trifocal_line_preSDP_mappingMatrix()\n','wt');
        fprintf(fid,'M=[\n');
        for it=1:Nt
            s=char(t(it));
            tok1=regexp(s,'R([0-9]+)(\*|\^)','tokens');
            if isempty(tok1)
                %this is the position of constant coefficients
                m=it; 
            else
                idx1=str2double(tok1{1}{1});
                tok2=regexp(s,'\*R([0-9]+)','tokens');
                if isempty(tok2)
                    %this is the position of a squared coefficient
                    idx2=idx1;
                else
                    idx2=str2double(tok2{1}{1});
                end
                M(it,sub2ind([18 18],idx1,idx2))=1;
                M(it,sub2ind([18 18],idx2,idx1))=1;
            end
            fprintf(fid,'\t');
            fprintf(fid,'%d ',M(it,:));
            fprintf(fid,'\n');
        end
        fprintf(fid,'];\n');
        fprintf(fid,'m=%d;\n',m);
        fclose(fid);
end
keyboard

function [Q,R]=flipqr(A)
[Q,R]=qr(A);
Q=fliplr(Q);
R=flipud(R);

function checkConstraint(R12,R23,data)
NLines=size(data.l3d,3);
for iLine=1:NLines
    l1=data.l2d(:,iLine,1);
    l2=data.l2d(:,iLine,2);
    l3=data.l2d(:,iLine,3);
    disp(l1'*R12*hat(l2)*R23*l3)
end
