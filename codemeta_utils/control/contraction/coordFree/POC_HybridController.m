% Perform various simulation runs using hybrid system
function [] = POC_HybridController(filename)

% Load/Define parameters
beta = [];
if ~isempty(filename) && exist(filename,'file')
    load(filename);
else
    load('data_ACC2019.mat');
end
J = diag([5;2;1]);
TFinal = 6;
w_rand = cnormalize(randn(3,1)); %random tangent vector

% Tests
R0(:,:,1) = rot_expVec(eye(3),pi*w_rand); % pi dist away from the identity
w0(:,1) = 10*cnormalize(randn(3,1));
R0(:,:,2) = rot_expVec(eye(3),pi/2*w_rand); % pi/2 dist away from the identity
w0(:,2) = 1/2*(sqrt(1/2*abs(kd)*maxD^2)+1)*w_rand;
R0(:,:,3) = rot_expVec(eye(3),maxD*w_rand); % pi dist away from the identity
w0(:,3) = magW*w_rand;
R0(:,:,4) = rot_expVec(eye(3),pi*w_rand); % pi dist away from the identity
w0(:,4) = zeros(3,1);
R0(:,:,5) = eye(3);
w0(:,5) = sqrt(1/2*abs(kd)*(pi-.1)^2)*w_rand;
R0(:,:,6) = rot_expVec(eye(3),(pi-1e-2)*w_rand); % pi dist away from the identity
w0(:,6) = 1*w_rand;
testName = ["Max Distance, Large Velocity","Inside Contraction Region", ...
    "Leaving Contraction Region","At pi, no vel", ...
    "At Identitiy with max speed for V","Crossing pi, low speed"];
for i = 1:1%size(w0,2)
    figure
    [R,w,t,u,state, RRef, wRef] = rot_SimulateDyn(R0(:,:,i),w0(:,i),TFinal, kd, kv, J, beta,...
        [m_contract(1,1);m_contract(1,2);m_contract(2,2)]);
    title(testName(i));
    [V, Vdot] = rotBundle_Lyapunov(R,w,kd,kv);
    figure
    plotFigures(R,w,t,u,kd,kv,V,Vdot,m_contract,beta,state,magW, RRef, wRef);
    sgtitle(testName(i))
end
end

function [] = plotFigures(R,w,t,u,kd,kv,V,Vdot,m_contract,beta,state,magW, RRef, wRef)
    n = 2; nn = 3; % define the matrix of figures
    % Plot States/Controller
    ax1 = subplot(n,nn,1);
    plotc(t,state,'o','LineWidth',3);
    title('State/Controller')
    % Plot Lyapunov to desired trajectory
    ax2 = subplot(n,nn,2);
    plotc(t,V,'o','LineWidth',3);
    title('Lyapunov wrt Identity');
    % Plot Lyapunov wrt Reference trajectory
    ax3 = subplot(n,nn,3);
    iState2Idx = find(state==2);
    if ~isempty(iState2Idx)
        for i = 1:length(iState2Idx)
            Vtest(i) = rotBundle_Lyapunov(R(:,:,iState2Idx(i)),...
                w(:,iState2Idx(i)),kd,kv,...
                'rreference',RRef(:,:,iState2Idx(i)),...
                'wreference',wRef(:,iState2Idx(i)));
        end
        plotc(t(iState2Idx),Vtest,'o','LineWidth',3);
    end
    title('Lyapunov wrt Reference Trajectory')
    % Plot Reference position error
    ax4 = subplot(n,nn,4);
    RefDist_Error = vecnorm(rot_vee(RRef,rot_log(RRef,R))).^2;
    plotc(t,RefDist_Error,'o','LineWidth',3);
    title('Ref Distance Error')
    % Plot Reference velocity error
    ax5 = subplot(n,nn,5);
    RefVel_Error = vecnorm(wRef-w);
    plotc(t,RefVel_Error,'o','LineWidth',3);
    title('Ref Velocity Error');
    % Plot control
    ax6 = subplot(n,nn,6);
    plotc(t,u,'o','LineWidth',3);
    title('Control')
    
%     linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
end

function h=plotc(X, varargin)
h=plot(X, varargin{:});
set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)));
if nargout==0,
    clear h
end
end
