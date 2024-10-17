function [Xe,Xp] = IcraEstimation(y,sigma_y,E,lens,sigma_lens,X0,sigma_X0,t)


% ladder points (epsilon = 1/L = step size of numerical integration, eqns 9-14 in paper)
L = 5; 



% X0:  3 x nV, where N is the number of joints
% Sigma_X0: 3nV x 3nV

% noise parameters

% process noise for links
ssquaredang = 10^(-4); 

% initial uncertainty for velocity
ssquared_v0 = 10^(-4); 

% process noise for root
ssquaredvr = 10^(-2); 


%% parameters 

nE = size(E,1);
nV = nE+ 1; % # of joints



[~,idx] = sort(E(:,2));
E = E(idx,:);

% compute (directed) adjacency matrix
Adj = zeros(nV);
for ie  = 1:nE, Adj(E(ie,1), E(ie,2)) = 1; end


% matrix with parents
% parents of node i have 1 in 
% self edge included
Ap = eye(nV);
for ie = 1:nE, 
    Ap(E(ie,1), E(ie,2)) = 1;
    Ap(:, E(ie,2)) = Ap(:, E(ie,2)) | Ap(:, E(ie,1));
end

F = numel(t);% # of frames
Xe = zeros(3,nV,F); % returned 3D structure, updates
Xp = zeros(3,nV,F); % returned 3D structure, predictions
Xp(:,:,1) = X0;





% dimension of output
m = 3*nV ; 


%%  initialization of state vector

% state vector dimension
n = 7*(nV-1) + 3 ; 
% 3*nV positions ; 3*nE velocities; nE lengths 


x = zeros(n,1);
P = zeros(n);

% joint lengths
idxlens  = 6*nE+4:7*nE+3;
x(idxlens) = lens;
P(idxlens,idxlens) = sigma_lens;  


% positions from initial 3D estimate
x(1:3) = X0(:,1); % root position

for ie = 1:nE
    idxe  = 3*ie+1:3*ie+3;
    x(idxe) = X0(:,E(ie,2))-X0(:,E(ie,1));
    x(idxe) = x(idxe)/norm(x(idxe)); % unit vector
end

[V0,D0] = eig(sigma_X0);

P0 = zeros(size(sigma_X0)); % should be 3*nV x 3*nV
for iv = 1:size(V0,2) % cols of V0 should be equal to 3*nV
    v1 = V0(:,iv);
    for ie=1:nE
        idxe  = 3*ie+1:3*ie+3;
       v1(idxe) =  (eye(3)- (x(idxe)* x(idxe)'))*v1(idxe);
    end
    P0 = P0 + D0(iv)*(v1*v1');
end

P(1:3*nV,1:3*nV) = P0;

% velocities are initialized to 0 
% initialize velocity uncertainty
for ie = 1 : nE
    idxe  = 3*ie+1:3*ie+3;
    idxv  = 3*nV+3*ie-2:3*nV+3*ie;
    Vb =  null(x(idxe)');
    P(idxv,idxv) = ssquared_v0*( Vb(:,1)*Vb(:,1)' + Vb(:,2)*Vb(:,2)');
end


%% Main loop

for i=1:F

    i

   %% update step 
   
    ye = zeros(m,1); % h(hat(x))
    yt = zeros(m,1); % true measurement
    
    size(y)
    size(yt)
    
    yt  = y(:,i);

    

   % compute differential of output map
   C = zeros(m,n); 

   
   % get 3D positions from orientations and lengths
   Pos = zeros(3*nV,1);
   Pos(1:3) = x(1:3);
   
   for ie=1:nE
      Pos(3*E(ie,2)-2:3*E(ie,2)) =  Pos(3*E(ie,1)-2:3*E(ie,1)) + x(3+6*nE+ie)*x(3*E(ie,2)-2:3*E(ie,2));
   end
   % ATTN: assumptions
   % ATTN: E(ie,2) = ie + 1
   % ATTN: E(ie,1) < E(ie,2)
   
   
   
   
   % compute jacobian now
   
    ye(1:3) = Pos(1:3);
    C(1:3,1:3) =  eye(3);

    for iv = 2:nV
        
         ye(3*iv-2:3*iv) =  Pos(3*iv-2:3*iv); % estimated output
         
     
        C(3*iv-2:3*iv,1:3) =  eye(3);
        for jv =2:iv
            if Ap(jv,iv)
                    C(3*iv-2:3*iv, 3*jv-2:3*jv) = x(3+6*nE+jv-1) * eye(3);
                    C(3*iv-2:3*iv, 3+6*nE+jv-1) =   x(3*jv-2:3*jv);
            end
        end
        
        
    end

   
    
   % output noise covariance
    R = sigma_y; 

  % Measurement update
   Mn = P*C'/(C*P*C'+R); % Kalman gain
   tangent = Mn*(yt-ye); % dx, x_new = x + dx

   

   % update and covariance propagation
  P = (eye(size(P))-Mn*C)*P; 
  
  
  
  % propagate updated estimate and covariance   to the new tangent space
  [U,S,~] = svd(P);
 
    xnew = x + tangent; % non-Euclidean components will be fixed below
   
    
   
  for ie=1:nE
       idxp = 3*ie+1:3*ie+3;
       idxv =  3*nV+3*ie-2:3*nV+3*ie;
      
       % costly part from here
      [p,u,~,~,Ut] = ExponentialAndTransportTS(x(idxp),x(idxv),...
                                               tangent(idxp),tangent(idxv),...
                                             U([idxp idxv],:),L);
       % up to here                                   
           xnew(idxp)= p;
           xnew(idxv)= u;
           U([idxp idxv],:)= Ut;     
  end           
    

    
    P =  U*S*U';
  
    x = xnew;

    
    
      % get 3D positions from orientations and lengths
   Pos = zeros(3*nV,1);
   Pos(1:3) = x(1:3);
   
   for ie=1:nE
      Pos(3*E(ie,2)-2:3*E(ie,2)) =  Pos(3*E(ie,1)-2:3*E(ie,1)) + x(3+6*nE+ie)*x(3*E(ie,2)-2:3*E(ie,2));
   end
   Xe(:,:,i) = reshape(Pos,[3 nV]);
    

        %% model-based prediction 
    
        if  i<F

            % state  
            xnew  = x;

                for ie=1:nE
                     idxp = 3*ie+1:3*ie+3;
                     idxv =  3*nV+3*ie-2:3*nV+3*ie;
                     xnew(idxp)=  ExpMapSphere(x(idxp),x(idxv)  );
                    xnew(idxv)= ParallelTranspot(x(idxp),xnew(idxp),x(idxv)); 
                end        


                    % covariance  
                    Q = zeros(n); % process noise covariance
                    A = zeros(n); % jacobian of f(x_t,w_t)

                    A(1:3,1:3) = eye(3);
                    A(6*nE+4:n,6*nE+4:n) = eye(n-6*nE-4+1);

                 Q(1:3,1:3) = ssquaredvr*eye(3); % process noise for root
                 
                       for ie=1:nE
                        idxp = 3*ie+1:3*ie+3;
                        idxv =  3*nV+3*ie-2:3*nV+3*ie;

                            idxie = [idxp idxv];
                            A(idxie,idxie) =  DifferentialOfFsingle( x(idxp),x(idxv));      

                            V = null(xnew(idxp)'); % ATTN: null of prediction

                            Q(idxv,idxv)  = ssquaredang*( (V(:,1)*V(:,1)') + (V(:,2)*V(:,2)'));   % process noise

                       end

                       
                            P = A*P*A'+ Q ;     % P[n+1|n]

                x = xnew;
                
                
            % get 3D positions from orientations and lengths
            Pos = zeros(3*nV,1);
            Pos(1:3) = x(1:3);

            for ie=1:nE
              Pos(3*E(ie,2)-2:3*E(ie,2)) =  Pos(3*E(ie,1)-2:3*E(ie,1)) + x(3+6*nE+ie)*x(3*E(ie,2)-2:3*E(ie,2));
            end
            Xp(:,:,i+1) = reshape(Pos,[3 nV]);
                
                
                
                

        end     % if i < F 

    
        
end     % for i = : F