 clear;clc;
T=1000;                 %Time-serie length

x_true=zeros(2,T);      %Initialization of the vector of true values
x_true(:,1)=[0;0];        %Initialization of the vector of true values at t=N
y=zeros(1,T);           %Initialization of the vector of observations
Qr    = [0.1 0.125 0.15 0.20 0.25 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1 1.5 2 2.5 3 10 100];
sigma_w   = zeros(1,length(Qr));
s_sigma_w = zeros(1,length(Qr));
%for i = 1:length(Qr)
sW_AR = 1;               % standard deviation, sigma_(AR,w)
sW_LL = 1e-12;
sV = 0.01;
%sV    = sW_AR./Qr(i);
phi   = 0.8;
R     = sV^2;
rng('default')
N     = 10000;
A_LL  = 1;
A_AR  = phi;
A     = blkdiag(A_LL,A_AR);
Q_LL  = sW_LL^2;
Q_AR  = sW_AR^2;
Q     = blkdiag(Q_LL,Q_AR);

EX_S = zeros(1,T,N);
PX_S = zeros(1,T,N);
index_onestd   = cell(T,1);
index_twostd   = cell(T,1);
index_thirdstd = cell(T,1);
ratio_onestd   = zeros(1,T);
ratio_twostd   = zeros(1,T);
ratio_thirdstd = zeros(1,T);
sW_AR         = zeros(1,N);
%% Data simulation
 for n = 1:N
    % Truncated Gaussian
    sw2=2;
    mw2=5;
%     D = 1-0.5*(1+erf(-mw2/(sqrt(2)*sw2)));
%     N = (1/sqrt(2*pi))*exp(-mw2^2/(2*sw2^2));
%     mw2_new = mw2 + sw2*(N/D);
%     var_new = sw2^2*(1-((mw2*N)/(sw2*D))-(N/D)^2);
%     sw2_new = sqrt(var_new);
    % sigma initialised from prior knowledge
     w2hat_sample = normrnd(mw2,sw2);
     while w2hat_sample<0
        w2hat_sample = normrnd(mw2,sw2);
     end
      sW_AR(n) = sqrt(w2hat_sample);
      Q          = blkdiag(sW_LL^2, sW_AR(n)^2);        %sW_AR for non-initialisation
    %YT = zeros(T,1);

    for t = 2:T
        x_true(:,t) = mvnrnd(A*x_true(:,t-1),Q);
        YT(:,t)     = x_true(1,t) + x_true(2,t) + normrnd(0,sV);
    end
    %plot(YT)
    %% State estimation
    EX=zeros(4,T);                               %   X  = [X^LL  X^AR  W  hat(W^2)]
    EX(:,1)=[0 0 nan mw2]';
    PX(:,:,1)=diag([1,1,nan,sw2^2]);

    %% Prediction
    for t=2:T
        Ep       = [EX(1,t-1);phi * EX(2,t-1); 0];              % mu_t|t-1
%         s_w_sq   = EX(4,t-1)+ PX(4,4,t-1)./(4*EX(4,t-1))+ 0.5*PX(4,4,t-1)*(1/(16*EX(4,t-1)^3)) + PX(4,4,t-1)^2./(64*EX(4,t-1)^3) - 2*sqrt(EX(4,t-1))*(PX(4,4,t-1)/8)*EX(4,t-1)^(-1.5);
%         if t>=4
%             s_w_sq = s_w_sq * (t/(t-2));
%         end
        %s_w_sq   = 2*(EX(4,t-1)+ 2*EX(4,t-1)^2 + PX(4,4,t-1));
        m_w2hat     = EX(4,t-1);
        s_w2hat     = sqrt(PX(4,4,t-1));
%         D_n         = 1-0.5*(1+erf(-m_w2hat/(sqrt(2)*s_w2hat)));
%         N_n         = (1/sqrt(2*pi))*exp(-m_w2hat^2/(2*s_w2hat^2));
%         m_w2hat_new = m_w2hat + s_w2hat*(N_n/D_n);
%         P_w2hat_new = s_w2hat^2*(1-((m_w2hat*N_n)/(s_w2hat*D_n))-(N_n/D_n)^2);
        s_w_sq      = m_w2hat;
        %PX(4,4,t-1) = P_w2hat_new;
%         if t>=3
%             s_w_sq = s_w_sq * (t/(t-2));
%         end
        Sp       = [ PX(1,1,t-1)       0                          0  ;
                        0         phi^2 * PX(2,2,t-1)+s_w_sq    s_w_sq ;        %  V(X)_t|t-1
                        0           s_w_sq                      s_w_sq  ];
        C        = [1 1 0 ];
        SY       = C*Sp*C'+R;
        SYX      = Sp*C';
        my       = C*Ep;
         K       = SYX/SY;
    %% Ist Update step:
       EX1  = Ep+SYX/SY*(YT(:,t)-my);
       %PX1  = (eye(2)-K*C)*Sp*(eye(2)-K*C)'+K*R*K';
       PX1  = Sp-SYX/SY*SYX';

       EX(1:3,t)     = EX1;
       PX(1:3,1:3,t) = PX1;

    %% 2nd Update step:
       s_w2_sq  = 2*(EX(4,t-1))^2;
       m_w2     = EX(3,t)^2+PX(3,3,t);
       s_w2     = 2*PX(3,3,t)^2+4*PX(3,3,t)*EX(3,t)^2;

       my1      = EX(4,t-1);
       SY1      = PX(4,4,t-1) + s_w2_sq + s_w2;
       SYX1     = PX(4,4,t-1);
       K1       = SYX1/SY1;

   %% Smoother Equations
       E_W2_pos      = m_w2;
       E_W2_prior    = my1;
       C_W2_W2hat    = SYX1;
       P_W2_prior    = 3*PX(4,4,t-1) + s_w2_sq;
       P_W2_pos      = s_w2;
       J             = C_W2_W2hat/P_W2_prior;
       EX(end,t)     = EX(4,t-1)  + J*(E_W2_pos - E_W2_prior);
       PX(end,end,t) = PX(4,4,t-1) + J^2*P_W2_pos - C_W2_W2hat^2/P_W2_prior;
%       PX(end,end,t) = PX(4,4,t-1) + J*(P_W2_pos - P_W2_prior)*J';
%        EX(end,t)        = EX(4,t-1)+SYX1/SY1*(m_w2-my1);
%        PX(end,end,t)    = PX(4,4,t-1)-SYX1/SY1*SYX1';
%    PX(end,end,t) = (1 - K1)'*PX(3,3,t-1)*(1-K1) + K1*s_w2*K1';

     end
      EX_S(:,:,n) = EX(end,:);
      PX_S(:,:,n) = squeeze(PX(end,end,:))';
      U1(n,:) = squeeze(EX_S(:,:,n) + sqrt(PX_S(:,:,n)));
      L1(n,:) = squeeze(EX_S(:,:,n) - sqrt(PX_S(:,:,n)));
      U2(n,:) = squeeze(EX_S(:,:,n) + 2*sqrt(PX_S(:,:,n)));
      L2(n,:) = squeeze(EX_S(:,:,n) - 2*sqrt(PX_S(:,:,n)));
      U3(n,:) = squeeze(EX_S(:,:,n) + 3*sqrt(PX_S(:,:,n)));
      L3(n,:) = squeeze(EX_S(:,:,n) - 3*sqrt(PX_S(:,:,n)));
      Tr(n,:) = sW_AR(n)^2.*ones(1,1000);
%       qqplot(squeeze(EX_S(:,:,2050)))
%       figure;
%       subplot(3,1,1)
%       plot(U1(1,:),'b');hold on;plot(L1(1,:),'g');hold on;plot(squeeze(EX_S(:,:,1)),'k');hold on;plot(sW_AR(1).^2*ones(1,1000),'--r');legend('Upper','Lower','Expected','True')
%       subplot(3,1,2)
%       plot(U2(1,:),'b');hold on;plot(L2(1,:),'g');hold on;plot(squeeze(EX_S(:,:,1)),'k');hold on;plot(sW_AR(1).^2*ones(1,1000),'--r')
%       subplot(3,1,3)
%       plot(U3(1,:),'b');hold on;plot(L3(1,:),'g');hold on;plot(squeeze(EX_S(:,:,1)),'k');hold on;plot(sW_AR(1).^2*ones(1,1000),'--r')
%
  end
% U1 = squeeze(EX_S + sqrt(PX_S));
% L1 = squeeze(EX_S - sqrt(PX_S));
% U2 = squeeze(EX_S + 2*sqrt(PX_S));
% L2 = squeeze(EX_S - 2*sqrt(PX_S));
% U3 = squeeze(EX_S + 3*sqrt(PX_S));
% L3 = squeeze(EX_S - 3*sqrt(PX_S));
% Tr = sW_AR^2.*ones(T,N);
% index_onestd   = cell(T,1);
% index_twostd   = cell(T,1);
% index_thirdstd = cell(T,1);
% ratio_onestd   = zeros(1,T);
% ratio_twostd   = zeros(1,T);
% ratio_thirdstd = zeros(1,T);
for i = 1:1000
index_onestd{i,1} = find(Tr(:,i) >= L1(:,i) & Tr(:,i) <= U1(:,i));
ratio_onestd(i)   = 100*size(index_onestd{i},1)/10000;
index_twostd{i,1} = find(Tr(:,i) >= L2(:,i) & Tr(:,i) <= U2(:,i));
ratio_twostd(i)   = 100*size(index_twostd{i},1)/10000;
index_thirdstd{i,1} = find(Tr(:,i) >= L3(:,i) & Tr(:,i) <= U3(:,i));
ratio_thirdstd(i)   = 100*size(index_thirdstd{i},1)/10000;
end
figure;
plot(ratio_onestd,'g');hold on;plot(ratio_twostd,'r');hold on;plot(ratio_thirdstd,'b');hold on
plot(95*ones(1,1000),'--k');hold on; plot(68*ones(1,1000),'--k');
xlabel('t');
ylabel('C.I.');
legend('one std','two std','third std');
%% Mean values CI
mean(ratio_onestd)
mean(ratio_twostd)
mean(ratio_thirdstd)
estim_sigma_w   = EX(4,t);
estim_sigma_w_P = PX(4,4,t);
% sigma_w(i) = estim_sigma_w;
% s_sigma_w(i) = estim_sigma_w_P;
% %end
% %MSE = mean((sW - estim_sigma_w).^2);
%
%
% % ratio = [0.25 0.5 0.75 0.85 1 1.5 2 2.5 3 10 100];
% % sigma_w = [0.0026 0.0030 0.0035 0.0037 0.0039 0.0043 0.0045 0.0048 0.0053 0.0096 0.0114];
% % s_sigma_w = [1.63e-04 1.87e-04 2.2e-04 2.32e-04 2.45e-04 2.64e-04 2.68e-04 2.77e-04 2.92e-04 4.53e-04 6.21e-04];
% %  ratio     = [0.1       0.125      0.15     0.20     0.25           0.30     0.40    0.5     0.60    0.70    0.80    0.90    0.95         1      1.5          2       2.5       3         10    ];
% %  sigma_w   = [4.74     4.4684      4.1392   3.4507   2.8533        2.3887   1.7825   1.4475   1.26    1.1542  1.0944  1.0607  1.0499    1.0418  1.02         1.0202    1.02    1.02     1.015    ];
% %  s_sigma_w = [0.6056   0.5331      0.4525   0.3062   0.2032        0.1380   0.0716   0.0716   0.0295  0.0219  0.0171  0.0127  0.0127    0.0116  0.0061      0.0042     0.0034    0.0030     0.0021];
%
% figure;
% patch([Qr,fliplr(Qr)],[sigma_w+s_sigma_w,fliplr(sigma_w-s_sigma_w)],'g','FaceAlpha',0.2,'EdgeColor','none');hold on;
% plot(Qr,sigma_w,'k');hold on; plot(Qr,repmat(sW_AR^2,[1,length(Qr)]),'-.r','Linewidth',1.5);
% xlabel('Q/R ratio');
% ylabel('sigma^2_w');ylim([0 0.1]);
% set(gca, 'XScale', 'log');
%%  Plotting Sigma_v
t  = 1:length(EX);
xw = EX(4,t);
sX = permute(sqrt(PX(4,4,t)),[1,3,2]);
figure;
plot(t,repmat(sW_AR^2,[1,length(EX)]),'-.r','Linewidth',1.5)
hold on;
patch([t,fliplr(t)],[xw+sX,fliplr(xw-sX)],'g','FaceAlpha',0.2,'EdgeColor','none')
hold on;
plot(t,xw,'k')
hold off
ylim([sW_AR^2-sW_AR/2,sW_AR^2+sW_AR^2/2])

xlabel('$Epoch$','Interpreter','latex')
ylabel('$sv$','Interpreter','latex')
