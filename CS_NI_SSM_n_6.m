% Square-root NI
% Created by: Bhargob Deka
% Dimension 6
% [d1(1)   d2(2)    d3(3)    d4(4)  d5(5)   d6(6)];  [  c1(7)    c2(8)    c3(9)    c4(10)
% c5(11)   c6(12)  c7(13)   c8(14)  c9(15)  c10(16)    c11(17)  c12(18)  c13(19)   c14(20) c15(21)]
syms d1 d2 d3 d4 d5 d6 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 real
A = [d1 c1 c2   c3   c4     c5;...
     0  d2 c6   c7   c8     c9;...
     0   0 d3  c10   c11   c12;...
     0   0  0   d4   c13   c14;...
     0   0  0   0    d5    c15;...
     0   0  0   0    0     d6];
 Sigma_W = A*A';
 
%  Sigma_W =
%  
% [ c1^2 + c2^2 + c3^2 + c4^2 + c5^2 + d1^2, c2*c6 + c3*c7 + c4*c8 + c5*c9 + c1*d2, c3*c10 + c4*c11 + c5*c12 + c2*d3,    c4*c13 + c5*c14 + c3*d4,   c5*c15 + c4*d5,  c5*d6]
% [   c2*c6 + c3*c7 + c4*c8 + c5*c9 + c1*d2,      c6^2 + c7^2 + c8^2 + c9^2 + d2^2, c7*c10 + c8*c11 + c9*c12 + c6*d3,    c8*c13 + c9*c14 + c7*d4,   c9*c15 + c8*d5,  c9*d6]
% [        c3*c10 + c4*c11 + c5*c12 + c2*d3,      c7*c10 + c8*c11 + c9*c12 + c6*d3,     c10^2 + c11^2 + c12^2 + d3^2, c11*c13 + c12*c14 + c10*d4, c12*c15 + c11*d5, c12*d6]
% [                 c4*c13 + c5*c14 + c3*d4,               c8*c13 + c9*c14 + c7*d4,       c11*c13 + c12*c14 + c10*d4,       c13^2 + c14^2 + d4^2, c14*c15 + c13*d5, c14*d6]
% [                          c5*c15 + c4*d5,                        c9*c15 + c8*d5,                 c12*c15 + c11*d5,           c14*c15 + c13*d5,     c15^2 + d5^2, c15*d6]
% [                                   c5*d6,                                 c9*d6,                           c12*d6,                     c14*d6,           c15*d6,   d6^2]
%  

%% MNI
clear;clc;
T=1000;                 %Time-serie length
n_x     = 6;
n_w     = n_x;
n       = n_x*2;
n_w2hat = n_x*(n_x+1)/2;
x_true  = zeros(n_x,T);      %  Initialization of the vector of true values
y       = zeros(T,1);        %  Initialization of the vector of observations
sV      = 0.001;
R       = sV^2;

% rng('default')
%% A Matrix
A_LL  = 1;
Block_A = cell(1,n_x);
for i = 1:n_x
    Block_A{i} = A_LL;
end
A = blkdiag(Block_A{:});
%% Q matrix
%sW_LL   = [2.1 1.7 2.2 1.8 1.6 2.4];
sW_LL   = [0.7781    0.8882    0.7213    0.6455    0.2651    0.4418];
% sW_LL   = 2*ones(1,6);
Block_Q = cell(1,n_x);
for i = 1:n_x
    Block_Q{i} = sW_LL(i);
end
% sW_LL1 = 1; sW_LL2 = 1; 
% sW_LL3 = 1; sW_LL4 = 1;
% Q_LL    = {sW_LL1^2 sW_LL2^2 sW_LL3^2 sW_LL4^2};
Q       = blkdiag(Block_Q{:});
true_value = 0.7;
% % Off-Diagonal terms in Q
% sW_AR12 = 4.3187e-01;
% sW_AR13 = -1.8910e-01;
% sW_AR14 = -1.9650e-01;
% sW_AR15 = 3.9222e-02;
% sW_AR16 = -1.8258e-01;
% sW_AR23 = -3.7298e-01;
% sW_AR24 = 4.1932e-02; 
% sW_AR25 = 3.2295e-01;
% sW_AR26 = -4.9368e-02; 
% sW_AR34 = -3.3614e-01;
% sW_AR35 = 6.7629e-02;
% sW_AR36 = -2.1376e-01;
% sW_AR45 = -9.7507e-02;
% sW_AR46 = 3.0103e-01;
% sW_AR56 = 8.5939e-02;

% sW_AR12 = 0.5;
% sW_AR13 = 0.3;
% sW_AR14 = 0.2;
% sW_AR15 = 0.4;
% sW_AR16 = 0.2;
% sW_AR23 = 0.3;
% sW_AR24 = 0.2; 
% sW_AR25 = 0.3;
% sW_AR26 = 0.4; 
% sW_AR34 = 0.3;
% sW_AR35 = 0.2;
% sW_AR36 = 0.3;
% sW_AR45 = 0.2;
% sW_AR46 = 0.3;
% sW_AR56 = 0.2;

sW_AR12 = -0.1381;
sW_AR13 = 0.336;
sW_AR14 = -0.1325;
sW_AR15 = 0.2268;
sW_AR16 = 0.3762;
sW_AR23 = 0.1830;
sW_AR24 = -0.3719; 
sW_AR25 = 0.1079;
sW_AR26 = -0.1186; 
sW_AR34 = -0.2838;
sW_AR35 = 0.3304;
sW_AR36 = 0.33325;
sW_AR45 = 0.0151;
sW_AR46 = -0.1858;
sW_AR56 = 0.0909;
sW_cov  = [[sW_AR12 sW_AR13 sW_AR14 sW_AR15 sW_AR16]'; [sW_AR23 sW_AR24 sW_AR25 sW_AR26]'; [sW_AR34 sW_AR35 sW_AR36]'; [sW_AR45 sW_AR46]'; sW_AR56 ];

Q(1,2:6) = [sW_AR12 sW_AR13 sW_AR14 sW_AR15 sW_AR16];
Q(2:6,1) = [sW_AR12 sW_AR13 sW_AR14 sW_AR15 sW_AR16]';
Q(2,3:6) = [sW_AR23 sW_AR24 sW_AR25 sW_AR26];
Q(3:6,2) = [sW_AR23 sW_AR24 sW_AR25 sW_AR26]';
Q(3,4:6) = [sW_AR34 sW_AR35 sW_AR36];
Q(4:6,3) = [sW_AR34 sW_AR35 sW_AR36]';
Q(4,5:6) = [sW_AR45 sW_AR46];
Q(5:6,4) = [sW_AR45 sW_AR46]';
Q(5,6)   = sW_AR56;
Q(6,5)   = sW_AR56;

w       = chol(Q)'*randn(n_x,T);
R       = sV^2.*eye(n_x);
v       = chol(R)'*randn(n_x,T);
%rng(1,'philox')
% rand_seed=4;
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',rand_seed));  %Initialize random stream number based on clock

%% Data
% Real Data 
%load('NI_SSM_Case3_Testset01');
YT = zeros(T,n_x);
% Observation matrix    
C  = eye(n_x);
% x1_true(:,1)=x_true(1:2,1);
% x2_true(:,1)=x_true(3:4,1);
for t = 2:T
    %x_true(:,t) = mvnrnd(A*x_true(:,t-1),Q);
    x_true(:,t)  = A*x_true(:,t-1) + w(:,t-1);
    YT(t,:)      = C*x_true(:,t)   + v(:,t-1);
%     YT(t,2)      = x_true(2,t) + normrnd(0,sV);
%     YT(t,3)      = x_true(3,t) + normrnd(0,sV);
%     YT(t,4)      = x_true(4,t) + normrnd(0,sV);
end
% data.values(:,1)     = YT(:,1);
% data.values(:,2)     = YT(:,2);
% GT                   = load('Data_D_T_RealData_Test2.mat');
% data.timestamps      = GT.timestamps(1:1095);
% data.labels          = {'T01','T02'};
% save('Data_T1_T2_NI_v1.mat',  '-struct', 'data');
% %YT(1305+1:end)       = nan;
% plot(YT(:,1),'k');hold on; plot(YT(:,2),'r')
%% State estimation
E_Lw        = zeros(21,T);
P_Lw        = zeros(21,21,T);
EL          = [0.5*ones(6,1);zeros(15,1)];
% PL          = 0.15*ones(21,1);
PL          = [0.1*ones(6,1); 0.05*ones(15,1)];
E_Lw(:,1)   = EL;                 
P_Lw(:,:,1) = diag(PL);

EX          = zeros(n_x+n_w+n_w2hat,T);
EX(:,1)     = [zeros(1,n_x) nan*zeros(1,n_x) EL']';                      
PX(:,:,1)   = diag([0.1*ones(1,n_x),nan*zeros(1,n_x),PL']);               

%my   = zeros(T,2);
% SY  = zeros(T,2);

%% Prediction 
for t=2:T
    Ep                = [A*EX(1:n_x,t-1);zeros(n_x,1)];            % mu_t|t-1
    E_L               = E_Lw(:,t-1);
    P_L               = diag(P_Lw(:,:,t-1));
    [E_P,P_P,Cov_L_P] = Chol_scaled(E_L,P_L);
    s_w_sq            = E_P;
    n_V               = n_x;    n_Cov = n_w2hat-n_x;               % no of variance terms and covariance terms for W2hat
    s_wii             = s_w_sq(1:n_V);                             % variance terms 
    s_wij             = s_w_sq(n_V+1:end);                         % covariance terms
    Sx                = A*PX(1:n_x,1:n_x,t-1)*A';         
%% Adding the var(W) to the Hidden states:
for i = 1 :n_x
    Sx(i,i) = Sx(i,i) + s_wii(i);
end
%% Adding the covariances between the hidden states
i = 1;
j = 1;
k = 1;
while i <= n_x - 1
Sx(i,j+1) = Sx(i,j+1) + s_wij(k);
 j = j+1;
 k = k+1;
  if j == n_x
     i = i+1;
     j = i;
  end
end
Sx = triu(Sx)+triu(Sx,1)';
Sx = (Sx + Sx')/2;
%% Matrices for the Variance of W
% Variance of W
P_W = diag(s_wii);
i = 1; j = 1; k = 1;
% Fill the covariance matrix with the elements from the vector of
% covariance
while i <= n_x-1
    P_W(i,j+1) = s_wij(k);
    j = j+1;
    k = k+1;
    if j == n_x
        i = i+1;
        j = i;
    end
end
% Filling the lower triangular matrix using the upper
P_W = triu(P_W)+triu(P_W,1)';

% Covariance matrix at t|t-1
Sp  = [Sx P_W;P_W P_W];
    if any(isnan(Sp))
        check = 1;
    end
% Observation matrix    
C  = [eye(n_x) zeros(n_x)];

SY           = C*Sp*C'+R.*eye(n_x);   %BD
SYX          = Sp*C';
my           = C*Ep;
K            = SYX/(SY + 1e-08);
e            = YT(t,:)'-my;
    if isnan(YT(t))
        EX(:,t)   = [Ep;s_w_sq'];
        PX(:,:,t) = blkdiag(Sp,PX(end-2:end,end-2:end,t-1));
    else
        %% Ist Update step:
        EX1  = Ep+K*e;
        %PX1  = (eye(n)-K*C)*pinv(Sp)*(eye(n)-K*C)' + K*R*K';
        PX1  = Sp-K*SYX';
        %PX1  = (PX1+PX1')/2;
        
        EX(1:n,t)     = EX1;    % n = n_x*2 i.e., no of x + no of w
        PX(1:n,1:n,t) = PX1;
        
        % Collecting W|y
        EX_wy   = EX1(end-n_x+1:end,1);
        PX_wy   = PX1(end-n_x+1:end,end-n_x+1:end,1);
        m       = triu(true(size(PX_wy)),1);  % Finding the upper triangular elements
        cwiwjy  = PX_wy(m)';                  % covariance elements between W's after update
        P_wy    = diag(PX_wy);
        i = 1; j = 1; k = 1;
        while i <= n_x-1
        PX_wy(i,j+1) = cwiwjy(k);
        j = j+1;
        k = k+1;
            if j == n_x
                i = i+1;
                j = i;
            end
        end
        PX_wy  = triu(PX_wy)+triu(PX_wy,1)';   % Updated covariance matrix of W
        
        %% Indices
        % Creating the cells that will hold the covariance terms
        % cov(w_iw_j,w_kw_l) and the indices for each row 
        n_w2      = n_x*(n_x+1)/2;
        cov_wijkl = cell(1,n_w2-1);
        ind_wijkl = cell(1,n_w2-1);
        ind_cov   = cell(1,n_w2-1);
        i = 1; s = n_w2;
        while i <= n_w2-1
            cov_wijkl{i}       = zeros(1,s-1);
            cov_prior_wijkl{i} = zeros(1,s-1);
            ind_wijkl{i}       = repmat(zeros(1,s-1)',1,4);
            ind_cov{i}         = repmat(zeros(1,s-1)',1,4);
            ind_cov_prior{i}   = repmat(zeros(1,s-1)',1,4);
            i = i+1;
            s = s-1;
        end
        v = 1:1:n_x;
        n_w = n_x;
        I1 = [[[1:n_w]' [1:n_w]']; nchoosek(v,2)];
        I2 = I1;
        I1(end,:) = [];
        for i = 1:n_w2-1
        ind_wijkl{i}(:,1:2) = repmat(I1(i,:),size(ind_wijkl{i},1),1);
        end
        i = 1; j = 2;
        while i <= n_w2-1
            ind_wijkl{i}(:,3:4) = I2(j:end,:);
            j = j+1;
            i = i+1;
        end
        %% Getting the mean indices
        ind_mu = ind_wijkl;
        i = 1;
        while i <= n_w2-1
        ind_mu{i} = [ind_mu{i}(:,1) ind_mu{i}(:,4) ind_mu{i}(:,2) ind_mu{i}(:,3)];
        i = i+1;
        end
        %% Getting the covariance indices from the \Sigma_W matrix
        i = 1;
        while i <= n_w2-1
              m_ijkl = ind_wijkl{i};  % matrix of all ijkl
              r      = size(m_ijkl,1); % no of rows
              for j = 1:r
                  ijkl = m_ijkl(j,:);  % vector ijkl for row j
                  ind_cov{i}(j,1) = n_w*(ijkl(3)-1)+ijkl(1);
                  ind_cov{i}(j,2) = n_w*(ijkl(4)-1)+ijkl(2);
                  ind_cov{i}(j,3) = n_w*(ijkl(3)-1)+ijkl(2);
                  ind_cov{i}(j,4) = n_w*(ijkl(4)-1)+ijkl(1);
              end
             i = i+1; 
        end
        
        %% 2nd Update step:
        % Computing W^2|y
        m_wii_y       = EX_wy.^2+diag(PX_wy);
        s_wii_y       = 2.*diag(PX_wy).^2+4.*diag(PX_wy).*EX_wy.^2;
        m_wiwj_y      = zeros(1,n_w2hat-n_x);
        s_wiwj_y      = zeros(1,n_w2hat-n_x);
        i = 1; j = 1; k = 1;
        while i <= n_x-1
            m_wiwj_y(k) = EX_wy(i)*EX_wy(j+1) + cwiwjy(k);
            s_wiwj_y(k) = P_wy(i)*P_wy(j+1) + cwiwjy(k)^2 + 2*cwiwjy(k)*EX_wy(i)*EX_wy(j+1) + P_wy(i)*EX_wy(j+1)^2 + P_wy(j+1)*EX_wy(i)^2;
            j = j+1;
            k = k+1;
            if j == n_x
                i = i+1;
                j = i;
            end
        end
        % Computing the covariance matrix \Sigma_W^2|y:
        i = 1;
        while i <= n_w2-1
            ind_C = ind_cov{i};
            ind_M = ind_mu{i};
            for j = 1:size(ind_C,1)
                cov_wijkl{i}(j) = cov1234(ind_C(j,:),ind_M(j,:),PX_wy,EX_wy);
            end
            n1 = (n_w2 - 1) - size(cov_wijkl{i},2);
                if n1 > 0
                    add_zeros = zeros(1,n1);
                else
                    add_zeros = [];
                end
            cov_wijkl{i} = [add_zeros cov_wijkl{i}];
            i = i+1;
        end
        %% Adding the variances and covariances of W^p to form \Sigma_W^p
        cov_wpy              = cell2mat(reshape(cov_wijkl,size(cov_wijkl,2),1));
        PX_wpy               = diag([s_wii_y' s_wiwj_y]);
        s_wpy                = zeros(size(PX_wpy,1));
        s_wpy(1:end-1,2:end) = cov_wpy;
        PX_wpy               = PX_wpy + s_wpy;
        PX_wpy               = triu(PX_wpy)+triu(PX_wpy,1)'; % adding the lower triangular matrix
        
        
        %% Creating E[W^p]
        EX_wpy = [m_wii_y' m_wiwj_y];

%% Computing prior mean and covariance matrix of Wp
        m_wsqhat    = E_P;
        s_wsqhat    = P_P;
        PX_wp       = zeros(size(s_wsqhat,1));
        ind12       = nchoosek(v,2);
        i = 1; j = 1;
        while i <= n_w2
            if i <= n_w
                PX_wp(i,i)   = 2*m_wsqhat(i)^2 + 3*s_wsqhat(i,i);
            else
                PX_wp(i,i)   = s_wsqhat(i,i) + (m_wsqhat(i)^2/(m_wsqhat(ind12(j,1))*m_wsqhat(ind12(j,2))+m_wsqhat(i)^2))*s_wsqhat(i,i)...
                    + m_wsqhat(ind12(j,1))*m_wsqhat(ind12(j,2)) + m_wsqhat(i)^2;
                j = j+1;
            end
            i = i+1;
        end
        %% Computing the indices for the covariances for prior \Sigma_W^2:
        i = 1;
        while i <= n_w2-1
              m_ijkl = ind_wijkl{i};  % matrix of all ijkl
              r      = size(m_ijkl,1); % no of rows
              v      = 1:1:n_x;
              for j = 1:r
                  ijkl = m_ijkl(j,:);  % vector ijkl for row j
                  if ijkl(1)==ijkl(3)
                      ind_cov_prior{i}(j,1) = find(v==ijkl(1));
                  else
                      diff = abs(ijkl(3) - ijkl(1));
                      ind_cov_prior{i}(j,1) = diff + n_x;
                  end
                  if ijkl(2)==ijkl(4)
                      ind_cov_prior{i}(j,2) = find(v==ijkl(2));
                  else
                      diff = abs(ijkl(4) - ijkl(2));
                      ind_cov_prior{i}(j,2) = diff + n_x;
                  end
                  if ijkl(1)==ijkl(4)
                      ind_cov_prior{i}(j,3) = find(v==ijkl(1));
                  else
                      diff = abs(ijkl(4) - ijkl(1));
                      ind_cov_prior{i}(j,3) = diff + n_x;
                  end
                  if ijkl(2)==ijkl(3)
                      ind_cov_prior{i}(j,4) = find(v==ijkl(2));
                  else
                      diff = abs(ijkl(3) - ijkl(2));
                      ind_cov_prior{i}(j,4) = diff + n_x;
                  end
                  
              end
             i = i+1; 
        end
        %% Computing the prior covariance matrix \Sigma_W^p
        i = 1;
        while i <= n_w2-1
            ind_Mpr = ind_cov_prior{i};
            for j = 1:size(ind_Mpr,1)
                cov_prior_wijkl{i}(j) = m_wsqhat(ind_Mpr(j,1)).*m_wsqhat(ind_Mpr(j,2)) + m_wsqhat(ind_Mpr(j,3)).*m_wsqhat(ind_Mpr(j,4));
            end
            n1 = (n_w2 - 1) - size(cov_prior_wijkl{i},2);
                if n1 > 0
                    add_zeros = zeros(1,n1);
                else
                    add_zeros = [];
                end
            cov_prior_wijkl{i} = [add_zeros cov_prior_wijkl{i}];
            i = i+1;
        end
        %% Adding the variances and covariances for Prior \Sigma_W^p
        cov_wp              = cell2mat(reshape(cov_prior_wijkl,size(cov_prior_wijkl,2),1));
        s_wp                = zeros(size(PX_wp,1));
        s_wp(1:end-1,2:end) = cov_wp;
        PX_wp               = PX_wp + s_wp;
        PX_wp               = triu(PX_wp)+triu(PX_wp,1)'; % adding the lower triangular matrix
        
        %% Creating Prior E[W^p]
        E_wp        = m_wsqhat;

        %% Smoother Equations
        E_Wp_prior    = E_wp;
        C_Wp_W2hat    = s_wsqhat;
        P_Wp_prior    = PX_wp;
        E_Wp_pos      = EX_wpy;
        P_Wp_pos      = PX_wpy;
        J             = C_Wp_W2hat/(P_Wp_prior+1e-08);
%         if ~all(all(J>1))
%             index    = find(J>1 & J < -1);
%             if ~isempty(index)
%                 J(index) = 1;
%             end
%         end
        EX(end-n_w2+1:end,t)                 = E_P   +      J*(E_Wp_pos' - E_Wp_prior);
        PX(end-n_w2+1:end,end-n_w2+1:end,t)  = P_P   +      J*(P_Wp_pos - P_Wp_prior)*J';
        E_Pw_y = EX(end-n_w2+1:end,t);
        P_Pw_y = PX(end-n_w2+1:end,end-n_w2+1:end,t);
        %% Converting from Pw to Lw
        Jc = Cov_L_P/(P_P+1e-08);
        E_Lw(:,t)   = E_Lw(:,t-1) + Jc*(E_Pw_y - E_P);
        P_Lw(:,:,t) = P_Lw(:,:,t-1) + Jc*(P_Pw_y - P_P)*Jc';
        
    end
end
[E_W_T, V_W_T, cov_W] = original_scale(E_Lw(:,T),P_Lw(:,:,T));

%% Plotting
%  Plotting Variances
E_W = zeros(21,length(EX));
V_W = zeros(21,length(EX));
for t = 1:length(EX)
    E_L                   = E_Lw(:,t);
    P_L                   = diag(P_Lw(:,:,t));
    [E_W(:,t),P_P(:,:,t)] = Chol_scaled(E_L,P_L);
    V_W(:,t)              = diag(P_P(:,:,t));
end
t  = 1:length(EX);
figure;
for i=1:n_x
    subplot(3,2,i)
    
    xw = E_W(i,t);
    sX = sqrt(V_W(i,t));
    plot(t,repmat(sW_LL(i),[1,length(EX)]),'-.r','Linewidth',1)
    hold on;
    patch([t,fliplr(t)],[xw+sX,fliplr(xw-sX)],'g','FaceAlpha',0.2,'EdgeColor','none')
    hold on;
    plot(t,xw,'k')
    hold off
    xlabel('$t$','Interpreter','latex')
    ylabel(['$\sigma^{\mathtt{AR}}' num2str(i) '$'],'Interpreter','latex')
    ylim([0,2.5])
end

% Plotting Covariances
figure;
for i=1:n_w2-n_w
    subplot(5,3,i)
    xw = E_W(n_w+i,t);
    sX = sqrt(V_W(n_w+i,t));
    plot(t,repmat(sW_cov(i,1),[1,length(EX)]),'-.r','Linewidth',1)
    hold on;
    patch([t,fliplr(t)],[xw+sX,fliplr(xw-sX)],'g','FaceAlpha',0.2,'EdgeColor','none')
    hold on;
    plot(t,xw,'k')
    hold off
    xlabel('$t$','Interpreter','latex')
    ylabel(['$\sigma^{\mathtt{AR}}' num2str(i) '$'],'Interpreter','latex')
    ylim([-1,1])
end
%% Functions
function V = cov1234(ind_C,ind_M,PX_wy,EX_wy)
            V = 2*PX_wy(ind_C(1)).*PX_wy(ind_C(2)) + 2*PX_wy(ind_C(3)).*EX_wy(ind_M(1)).*EX_wy(ind_M(2)) + 2*PX_wy(ind_C(4)).*EX_wy(ind_M(3)).*EX_wy(ind_M(4));
end

function [E_P, P_P, Cov_L_P] = Chol_scaled(E_L,P_L)
%E_L is 21x1 vector of expected values
%P_L is 21x1 vector of variance terms
% E_L = [ones(6,1);zeros(15,1)];
% P_L = ones(21,1);
P_P     = zeros(21,1);
E_P     = zeros(21,1);
Cov_L_P = zeros(21,21);

E_P(1)  =   Mean_GMA(E_L,P_L,[7,7,7])+Mean_GMA(E_L,P_L,[8,8,8])+Mean_GMA(E_L,P_L,[9,9,9])+Mean_GMA(E_L,P_L,[10,10,10])+Mean_GMA(E_L,P_L,[11,11,11])+Mean_GMA(E_L,P_L,[1,1,1]);
E_P(2)  =   Mean_GMA(E_L,P_L,[12,12,12])+Mean_GMA(E_L,P_L,[13,13,13])+Mean_GMA(E_L,P_L,[14,14,14])+Mean_GMA(E_L,P_L,[15,15,15])+Mean_GMA(E_L,P_L,[2,2,2]);
E_P(3)  =   Mean_GMA(E_L,P_L,[16,16,16])+Mean_GMA(E_L,P_L,[17,17,17])+Mean_GMA(E_L,P_L,[18,18,18])+Mean_GMA(E_L,P_L,[3,3,3]);
E_P(4)  =   Mean_GMA(E_L,P_L,[19,19,19])+Mean_GMA(E_L,P_L,[20,20,20])+Mean_GMA(E_L,P_L,[4,4,4]);
E_P(5)  =   Mean_GMA(E_L,P_L,[21,21,21]) + Mean_GMA(E_L,P_L,[5,5,5]);
E_P(6)  =   Mean_GMA(E_L,P_L,[6,6,6]);

E_P(7)  =   Mean_GMA(E_L,P_L,[8,12,[]]) + Mean_GMA(E_L,P_L,[9,13,[]]) + Mean_GMA(E_L,P_L,[10,14,[]]) + Mean_GMA(E_L,P_L,[11,15,[]]) + Mean_GMA(E_L,P_L,[7,2,[]]);
E_P(8)  =   Mean_GMA(E_L,P_L,[9,16,[]]) + Mean_GMA(E_L,P_L,[10,17,[]]) + Mean_GMA(E_L,P_L,[11,18,[]]) + Mean_GMA(E_L,P_L,[8,3,[]]);
E_P(9)  =   Mean_GMA(E_L,P_L,[10,19,[]]) + Mean_GMA(E_L,P_L,[11,20,[]]) + Mean_GMA(E_L,P_L,[9,4,[]]) ;
E_P(10)  =  Mean_GMA(E_L,P_L,[11,21,[]]) + Mean_GMA(E_L,P_L,[10,5,[]]);
E_P(11)  =   Mean_GMA(E_L,P_L,[11,6,[]]);
% % [d1(1)   d2(2)    d3(3)    d4(4)  d5(5)   d6(6)];  [  c1(7)    c2(8)    c3(9)    c4(10)
% c5(11)   c6(12)  c7(13)   c8(14)  c9(15)  c10(16)    c11(17)  c12(18)  c13(19)   c14(20) c15(21)]

% c7*c10 + c8*c11 + c9*c12 + c6*d3,    c8*c13 + c9*c14 + c7*d4,  c9*c15+c8*d5,  c9*d6
E_P(12)  =  Mean_GMA(E_L,P_L,[13,16,[]])+ Mean_GMA(E_L,P_L,[14,17,[]])  + Mean_GMA(E_L,P_L,[15,18,[]]) + Mean_GMA(E_L,P_L,[12,3,[]]);
E_P(13)  =  Mean_GMA(E_L,P_L,[14,19,[]])+ Mean_GMA(E_L,P_L,[15,20,[]])  + Mean_GMA(E_L,P_L,[13,4,[]]);
E_P(14)  =  Mean_GMA(E_L,P_L,[15,21,[]])+ Mean_GMA(E_L,P_L,[14,5,[]]);
E_P(15)  =  Mean_GMA(E_L,P_L,[15,6,[]]);
% c11*c13 + c12*c14 + c10*d4, c12*c15+c11*d5, c12*d6
E_P(16)  =  Mean_GMA(E_L,P_L,[17,19,[]])+ Mean_GMA(E_L,P_L,[18,20,[]])  + Mean_GMA(E_L,P_L,[16,4,[]]);
E_P(17)  =  Mean_GMA(E_L,P_L,[18,21,[]])+ Mean_GMA(E_L,P_L,[17,5,[]]);
E_P(18)  =  Mean_GMA(E_L,P_L,[18,6,[]]);
% c14*c15+c13*d5, c14*d6
E_P(19)  =  Mean_GMA(E_L,P_L,[20,21,[]])+ Mean_GMA(E_L,P_L,[19,5,[]]);
E_P(20)  =  Mean_GMA(E_L,P_L,[20,6,[]]);
% c15*d6
E_P(21)  =  Mean_GMA(E_L,P_L,[21,6,[]]);
%% Variance
% % [d1(1)   d2(2)    d3(3)    d4(4)  d5(5)   d6(6)];  [  c1(7)    c2(8)    c3(9)    c4(10)
% c5(11)   c6(12)  c7(13)   c8(14)  c9(15)  c10(16)    c11(17)  c12(18)  c13(19)   c14(20) c15(21)]

% c1^2 + c2^2 + c3^2 + c4^2 + c5^2 + d1^2
P_P(1,1)  =   Var12_GMA(E_L,P_L,[7,7])+Var12_GMA(E_L,P_L,[8,8])+Var12_GMA(E_L,P_L,[9,9])+Var12_GMA(E_L,P_L,[1,1]);
% c6^2 + c7^2 + c8^2 + c9^2 + d2^2
P_P(2,1)  =   Var12_GMA(E_L,P_L,[12,12])+Var12_GMA(E_L,P_L,[13,13])+Var12_GMA(E_L,P_L,[14,14])+Var12_GMA(E_L,P_L,[15,15])+Var12_GMA(E_L,P_L,[2,2]);
% c10^2 + c11^2 + c12^2 + d3^2
P_P(3,1)  =   Var12_GMA(E_L,P_L,[16,16])+ Var12_GMA(E_L,P_L,[17,17])+ Var12_GMA(E_L,P_L,[18,18])+ Var12_GMA(E_L,P_L,[3,3]);
% c13^2 + c14^2 + d4^2
P_P(4,1)  =   Var12_GMA(E_L,P_L,[19,19]) + Var12_GMA(E_L,P_L,[20,20]) + Var12_GMA(E_L,P_L,[4,4]);
% d5^2
P_P(5,1)  =   Var12_GMA(E_L,P_L,[21,21]) + Var12_GMA(E_L,P_L,[5,5]);
% d6^2
P_P(6,1)  =   Var12_GMA(E_L,P_L,[6,6]);

% c2*c6 + c3*c7 + c4*c8 + c5*c9 + c1*d2
P_P(7,1)  =   Var12_GMA(E_L,P_L,[8,12]) + Var12_GMA(E_L,P_L,[9,13]) + Var12_GMA(E_L,P_L,[10,14]) + Var12_GMA(E_L,P_L,[11,15]) + Var12_GMA(E_L,P_L,[7,2]);
% c3*c10 + c4*c11 + c5*c12 + c2*d3
P_P(8,1)  =   Var12_GMA(E_L,P_L,[9,16]) + Var12_GMA(E_L,P_L,[10,17]) + Var12_GMA(E_L,P_L,[11,18]) + Var12_GMA(E_L,P_L,[8,3]);
% c4*c13 + c5*c14 + c3*d4
P_P(9,1)  =   Var12_GMA(E_L,P_L,[10,19]) + Var12_GMA(E_L,P_L,[11,20]) + Var12_GMA(E_L,P_L,[9,4]) ;
% c5*c15+c4*d5
P_P(10,1) =   Var12_GMA(E_L,P_L,[11,21]) + Var12_GMA(E_L,P_L,[10,5]);
% c5*d6
P_P(11,1) =   Var12_GMA(E_L,P_L,[11,6]);

% c7*c10 + c8*c11 + c9*c12 + c6*d3,    c8*c13 + c9*c14 + c7*d4,  c9*c15+c8*d5,  c9*d6
P_P(12,1)  =  Var12_GMA(E_L,P_L,[13,16])+ Var12_GMA(E_L,P_L,[14,17])  + Var12_GMA(E_L,P_L,[15,18]) + Var12_GMA(E_L,P_L,[12,3]);
P_P(13,1)  =  Var12_GMA(E_L,P_L,[14,19])+ Var12_GMA(E_L,P_L,[15,20])  + Var12_GMA(E_L,P_L,[13,4]);
P_P(14,1)  =  Var12_GMA(E_L,P_L,[15,21])+ Var12_GMA(E_L,P_L,[14,5]);
P_P(15,1)  =  Var12_GMA(E_L,P_L,[15,6]);

% c11*c13 + c12*c14 + c10*d4, c12*c15+c11*d5, c12*d6
P_P(16,1)  =  Var12_GMA(E_L,P_L,[17,19])+ Var12_GMA(E_L,P_L,[18,20])  + Var12_GMA(E_L,P_L,[16,4]);
P_P(17,1)  =  Var12_GMA(E_L,P_L,[18,21])+ Var12_GMA(E_L,P_L,[17,5]);
P_P(18,1)  =  Var12_GMA(E_L,P_L,[18,6]);
% c14*c15+c13*d5, c14*d6
P_P(19,1)  =  Var12_GMA(E_L,P_L,[20,21])+ Var12_GMA(E_L,P_L,[19,5]);
P_P(20,1)  =  Var12_GMA(E_L,P_L,[20,6]);
% c15*d6
P_P(21,1)  =  Var12_GMA(E_L,P_L,[21,6]);
P_P        =  diag(P_P);
% Covariance
% [ c1^2 + c2^2 + c3^2 + c4^2 + c5^2 + d1^2, c2*c6 + c3*c7 + c4*c8 + c5*c9 + c1*d2, c3*c10 + c4*c11 + c5*c12 + c2*d3,    c4*c13 + c5*c14 + c3*d4,  c5*d5,  c5*d6]
% [   c2*c6 + c3*c7 + c4*c8 + c5*c9 + c1*d2,      c6^2 + c7^2 + c8^2 + c9^2 + d2^2, c7*c10 + c8*c11 + c9*c12 + c6*d3,    c8*c13 + c9*c14 + c7*d4,  c9*d5,  c9*d6]
% [        c3*c10 + c4*c11 + c5*c12 + c2*d3,      c7*c10 + c8*c11 + c9*c12 + c6*d3,     c10^2 + c11^2 + c12^2 + d3^2, c11*c13 + c12*c14 + c10*d4, c12*d5, c12*d6]
% [                 c4*c13 + c5*c14 + c3*d4,               c8*c13 + c9*c14 + c7*d4,       c11*c13 + c12*c14 + c10*d4,       c13^2 + c14^2 + d4^2, c14*d5, c14*d6]
% [                                   c5*d5,                                 c9*d5,                           c12*d5,                     c14*d5,   d5^2,  d5*d6]
% [                                   c5*d6,                                 c9*d6,                           c12*d6,                     c14*d6,  d5*d6,   d6^2]

 % Column 1 -- c1^2 + c2^2 + c3^2 + c4^2 + c5^2 + d1^2
Cov_L_P(1,1)      = Cov123_GMA(E_L,P_L,[1 1 1]);
Cov_L_P(2:6,1)    = zeros(5,1); 
Cov_L_P(7,1)      = Cov123_GMA(E_L,P_L,[7 7 7]);
Cov_L_P(8,1)      = Cov123_GMA(E_L,P_L,[8 8 8]);
Cov_L_P(9,1)      = Cov123_GMA(E_L,P_L,[9 9 9]);
Cov_L_P(10,1)     = Cov123_GMA(E_L,P_L,[10 10 10]);
Cov_L_P(11,1)     = Cov123_GMA(E_L,P_L,[11 11 11]);
Cov_L_P(12:21,1)  = zeros(10,1);

% Column 2 -- c6^2 + c7^2 + c8^2 + c9^2 + d2^2
Cov_L_P(1,2)      =  0;
Cov_L_P(2,2)      =  Cov123_GMA(E_L,P_L,[2 2 2]);
Cov_L_P(3:6,2)    =  zeros(4,1);
Cov_L_P(7:11,2)   =  zeros(5,1);
Cov_L_P(12,2)     =  Cov123_GMA(E_L,P_L,[12 12 12]);
Cov_L_P(13,2)     =  Cov123_GMA(E_L,P_L,[13 13 13]);
Cov_L_P(14,2)     =  Cov123_GMA(E_L,P_L,[14 14 14]);
Cov_L_P(15,2)     =  Cov123_GMA(E_L,P_L,[15 15 15]);
Cov_L_P(16:21,2)  =  zeros(6,1);

% Column 3 -- c10^2 + c11^2 + c12^2 + d3^2
Cov_L_P(1,3)     = 0;
Cov_L_P(2,3)     = 0;
Cov_L_P(3,3)     = Cov123_GMA(E_L,P_L,[3 3 3]);
Cov_L_P(4:6,3 )  =  zeros(3,1);
Cov_L_P(7:15,3)  = zeros(9,1);
Cov_L_P(16,3)    = Cov123_GMA(E_L,P_L,[16 16 16]);
Cov_L_P(17,3)    = Cov123_GMA(E_L,P_L,[17 17 17]);
Cov_L_P(18,3)    = Cov123_GMA(E_L,P_L,[18 18 18]);
% Column 4 -- c13^2 + c14^2 + d4^2
Cov_L_P(1,4)    = 0;
Cov_L_P(2,4)    = 0;
Cov_L_P(3,4)    = 0;
Cov_L_P(4,4)    = Cov123_GMA(E_L,P_L,[4 4 4]);
Cov_L_P(5:6,4)  = zeros(2,1);
Cov_L_P(7:18,4) = zeros(18-7+1,1);
Cov_L_P(19,4)   = Cov123_GMA(E_L,P_L,[19 19 19]);
Cov_L_P(20,4)   = Cov123_GMA(E_L,P_L,[20 20 20]);

% Column 5 -- d5^2
Cov_L_P(1:4,5)  = zeros(4,1);
Cov_L_P(5,5)    = Cov123_GMA(E_L,P_L,[5 5 5]);
Cov_L_P(6,5)    = 0;
Cov_L_P(7:21,5) = zeros(21-7+1,1);
% Column 6 -- d6^2
Cov_L_P(1:5,6)  = zeros(5,1);
Cov_L_P(6,6)    = Cov123_GMA(E_L,P_L,[6 6 6]);
Cov_L_P(7:21,6) = zeros(21-7+1,1);
% Column 7 -- c2*c6 + c3*c7 + c4*c8 + c5*c9 + c1*d2
Cov_L_P(1,7)   = 0;
Cov_L_P(2,7)   = Cov123_GMA(E_L,P_L,[2 7 2]);
Cov_L_P(3:6,7) = zeros(4,1);

Cov_L_P(7,7)   = Cov123_GMA(E_L,P_L,[7 7 2]);
Cov_L_P(8,7)   = Cov123_GMA(E_L,P_L,[8 8 12]);
Cov_L_P(9,7)   = Cov123_GMA(E_L,P_L,[9 9 13]);
Cov_L_P(10,7)  = Cov123_GMA(E_L,P_L,[10 10 14]);
Cov_L_P(11,7)  = Cov123_GMA(E_L,P_L,[11 11 15]);
Cov_L_P(12,7)  = Cov123_GMA(E_L,P_L,[12 8 12]);
Cov_L_P(13,7)  = Cov123_GMA(E_L,P_L,[13 9 13]);
Cov_L_P(14,7)  = Cov123_GMA(E_L,P_L,[14 10 14]);
Cov_L_P(15,7)  = Cov123_GMA(E_L,P_L,[15 11 15]);
Cov_L_P(16:21,7)  = zeros(21-16+1,1);

% Column 8 --  c3*c10 + c4*c11 + c5*c12 + c2*d3  
Cov_L_P(1,8)    = 0;
Cov_L_P(2,8)    = 0;
Cov_L_P(3,8)    = Cov123_GMA(E_L,P_L,[3 8 3]);
Cov_L_P(4:6,8)  = zeros(3,1);
Cov_L_P(7,8)    = 0;
Cov_L_P(8,8)    = Cov123_GMA(E_L,P_L,[8 8 3]);
Cov_L_P(9,8)    = Cov123_GMA(E_L,P_L,[9 9 16]);
Cov_L_P(10,8)   = Cov123_GMA(E_L,P_L,[10 10 17]);
Cov_L_P(11,8)   = Cov123_GMA(E_L,P_L,[11 11 18]);
Cov_L_P(12:16,8)= zeros(16-12+1,1);
Cov_L_P(16,8)   = Cov123_GMA(E_L,P_L,[16 9 16]);
Cov_L_P(17,8)   = Cov123_GMA(E_L,P_L,[17 10 17]);
Cov_L_P(18,8)   = Cov123_GMA(E_L,P_L,[18 11 18]);
Cov_L_P(19:21,8)= zeros(21-19+1,1);

% Column 9 -- c4*c13 + c5*c14 + c3*d4
Cov_L_P(1,9)  = 0;
Cov_L_P(2,9)  = 0;
Cov_L_P(3,9)  = 0;
Cov_L_P(4,9)  = Cov123_GMA(E_L,P_L,[4 9 4]);
Cov_L_P(5,9)  = 0;
Cov_L_P(6,9)  = 0;
Cov_L_P(7,9)  = 0;
Cov_L_P(8,9)  = 0;
Cov_L_P(9,9)    = Cov123_GMA(E_L,P_L,[9 9 4]);
Cov_L_P(10,9)   = Cov123_GMA(E_L,P_L,[10 10 19]);
Cov_L_P(11,9)   = Cov123_GMA(E_L,P_L,[11 11 20]);
Cov_L_P(12:18,9)= zeros(18-12+1,1);
Cov_L_P(19,9)   = Cov123_GMA(E_L,P_L,[19 10 19]);
Cov_L_P(20,9)   = Cov123_GMA(E_L,P_L,[20 11 20]);
Cov_L_P(21,9)   = 0;
% Column 10 -- c5*c15+c4*d5
Cov_L_P(1:4,10)       = zeros(4,1);
Cov_L_P(5,10)         = Cov123_GMA(E_L,P_L,[5 10 5]);
Cov_L_P(6,10)         = 0;
Cov_L_P(7:9,10)       = zeros(3,1);
Cov_L_P(10,10)        = Cov123_GMA(E_L,P_L,[10 10 5]);
Cov_L_P(11,10)        = Cov123_GMA(E_L,P_L,[11 11 21]);
Cov_L_P(12:20,10)     = zeros(20-12+1,1);
Cov_L_P(21,10)        = Cov123_GMA(E_L,P_L,[21 11 21]);
% Column 11 -- c5*d6
Cov_L_P(1:5,11)      = zeros(5,1);
Cov_L_P(6,11)        = Cov123_GMA(E_L,P_L,[6 11 6]);
Cov_L_P(7:10,11)     = zeros(4,1);
Cov_L_P(11,11)       = Cov123_GMA(E_L,P_L,[11 11 6]);
Cov_L_P(12:21,11)    = zeros(21-12+1,1);

%% c7*c10 + c8*c11 + c9*c12 + c6*d3,    c8*c13 + c9*c14 + c7*d4,  c9*d5,  c9*d6
% Column 12 -- c7*c10 + c8*c11 + c9*c12 + c6*d3
Cov_L_P(1,12)     = 0;
Cov_L_P(2,12)     = 0;
Cov_L_P(3,12)     = Cov123_GMA(E_L,P_L,[3 12 3]);
Cov_L_P(4:6,12)   = zeros(3,1);
Cov_L_P(7:11,12)    = zeros(11-7+1,1);
Cov_L_P(12,12)    = Cov123_GMA(E_L,P_L,[12 12 3]);
Cov_L_P(13,12)    = Cov123_GMA(E_L,P_L,[13 13 16]);
Cov_L_P(14,12)    = Cov123_GMA(E_L,P_L,[14 14 17]);
Cov_L_P(15,12)    = Cov123_GMA(E_L,P_L,[15 15 18]);
Cov_L_P(16,12)    = Cov123_GMA(E_L,P_L,[16 13 16]);
Cov_L_P(17,12)    = Cov123_GMA(E_L,P_L,[17 14 17]);
Cov_L_P(18,12)    = Cov123_GMA(E_L,P_L,[18 15 18]);
Cov_L_P(19:21,12) = zeros(21-19+1,1);
% Column 13 -- c8*c13 + c9*c14 + c7*d4
Cov_L_P(1,13)     = 0;
Cov_L_P(2,13)     = 0;
Cov_L_P(3,13)     = 0;
Cov_L_P(4,13)     = Cov123_GMA(E_L,P_L,[4 13 4]);
Cov_L_P(5,13)     = 0;
Cov_L_P(6,13)     = 0;
Cov_L_P(7:12,13)  = zeros(12-7+1,1);
Cov_L_P(13,13)   = Cov123_GMA(E_L,P_L,[13 13 4]);
Cov_L_P(14,13)   = Cov123_GMA(E_L,P_L,[14 14 19]);
Cov_L_P(15,13)   = Cov123_GMA(E_L,P_L,[15 15 20]);
Cov_L_P(16:18,13)= zeros(18-16+1,1);
Cov_L_P(19,13)   = Cov123_GMA(E_L,P_L,[19 14 19]);
Cov_L_P(20,13)   = Cov123_GMA(E_L,P_L,[20 15 20]);
Cov_L_P(21,13)   = 0;
% Column 14 -- c9*c15+c8*d5
Cov_L_P(1:4,14)    = zeros(4,1);
Cov_L_P(5,14)      = Cov123_GMA(E_L,P_L,[5 14 5]);
Cov_L_P(6,14)      = 0;
Cov_L_P(7:13,14)   = zeros(13-7+1,1);
Cov_L_P(14,14)     = Cov123_GMA(E_L,P_L,[14 14 5]);
Cov_L_P(15,14)     = Cov123_GMA(E_L,P_L,[15 15 21]);
Cov_L_P(16:20,14)  = zeros(20-16+1,1);
Cov_L_P(21,14)     = Cov123_GMA(E_L,P_L,[21 15 21]);
% Column 15 -- c9*d6
Cov_L_P(1:5,15)      = zeros(5,1);
Cov_L_P(6,15)        = Cov123_GMA(E_L,P_L,[6 15 6]);
Cov_L_P(7:14,15)     = zeros(14-7+1,1);
Cov_L_P(15,15)       = Cov123_GMA(E_L,P_L,[15 15 6]);
Cov_L_P(16:21,15)    = zeros(21-16+1,1);
%% c11*c13 + c12*c14 + c10*d4, c12*d5, c12*d6
% Column 16 -- c11*c13 + c12*c14 + c10*d4
Cov_L_P(1,16)       = 0;
Cov_L_P(2,16)       = 0;
Cov_L_P(3,16)       = 0;
Cov_L_P(4,16)       = Cov123_GMA(E_L,P_L,[4 16 4]);
Cov_L_P(5,16)       = 0;
Cov_L_P(6,16)       = 0;
Cov_L_P(7:15,16)    = zeros(15-7+1,1);
Cov_L_P(16,16)      = Cov123_GMA(E_L,P_L,[16 16 4]);
Cov_L_P(17,16)      = Cov123_GMA(E_L,P_L,[17 17 19]);
Cov_L_P(18,16)      = Cov123_GMA(E_L,P_L,[18 18 20]);
Cov_L_P(19,16)      = Cov123_GMA(E_L,P_L,[19 17 19]);
Cov_L_P(20,16)      = Cov123_GMA(E_L,P_L,[20 18 20]);
Cov_L_P(21,16)      = 0;
% Column 17 -- c12*c15+c11*d5
Cov_L_P(1:4,17)     = zeros(4,1);
Cov_L_P(5,17)       = Cov123_GMA(E_L,P_L,[5 17 5]);
Cov_L_P(6,17)       = 0;
Cov_L_P(7:16,17)    = zeros(16-7+1,1);
Cov_L_P(17,17)      = Cov123_GMA(E_L,P_L,[17 17 5]);
Cov_L_P(18,17)      = Cov123_GMA(E_L,P_L,[18 18 21]);
Cov_L_P(19:20,17)   = zeros(2,1);
Cov_L_P(21,17)      = Cov123_GMA(E_L,P_L,[21 18 21]);
% Column 18 -- c12*d6
Cov_L_P(1:5,18)     = zeros(5,1);
Cov_L_P(6,18)       = Cov123_GMA(E_L,P_L,[6 18 6]);
Cov_L_P(7:17,18)    = zeros(17-7+1,1);
Cov_L_P(18,18)      = Cov123_GMA(E_L,P_L,[18 18 6]);
Cov_L_P(19:21,18)   = zeros(3,1);
% Column 19 -- c14*c15+c13*d5
Cov_L_P(1:4,19)     = zeros(4,1);
Cov_L_P(5,19)       = Cov123_GMA(E_L,P_L,[5 19 5]);
Cov_L_P(6,19)       = 0;
Cov_L_P(7:18,19)    = zeros(18-7+1,1);
Cov_L_P(19,19)      = Cov123_GMA(E_L,P_L,[19 19 5]);
Cov_L_P(20,19)      = Cov123_GMA(E_L,P_L,[20 20 21]);
Cov_L_P(21,19)      = Cov123_GMA(E_L,P_L,[21 20 21]);
% Column 20 -- c14*d6
Cov_L_P(1:5,20)     = zeros(5,1);
Cov_L_P(6,20)       = Cov123_GMA(E_L,P_L,[6 20 6]);
Cov_L_P(7:19,20)    = zeros(19-7+1,1);
Cov_L_P(20,20)      = Cov123_GMA(E_L,P_L,[20 20 6]);
Cov_L_P(21,20)      = 0;
% Column 21 -- c15*d6
Cov_L_P(1:5,21)     = zeros(5,1);
Cov_L_P(6,21)       = Cov123_GMA(E_L,P_L,[6 21 6]);
Cov_L_P(7:20,21)    = zeros(20-7+1,1);
Cov_L_P(21,21)      = Cov123_GMA(E_L,P_L,[21 21 6]);

end


function M = Mean_GMA(E_L,P_L,ind)
    if ind(1)==ind(2)
        M = E_L(ind(1))*E_L(ind(2))+P_L(ind(3));
    else
        M = E_L(ind(1))*E_L(ind(2));
    end
end
function V = Var12_GMA(E_L,P_L,ind)
    if ind(1)==ind(2)
        V = 2*P_L(ind(1))^2 + 4*P_L(ind(1))*E_L(ind(1))^2;
    else
        V = P_L(ind(1))*P_L(ind(2)) + P_L(ind(1))*E_L(ind(2))^2 + P_L(ind(2))*E_L(ind(1))^2;
    end
end
function C = Cov123_GMA(E_L,P_L,ind)
    if ind(1)==ind(2)==ind(3)
        C = 2*P_L(ind(1))*E_L(ind(1));
    elseif ind(1)==ind(2)
        C = P_L(ind(1))*E_L(ind(3));
    elseif ind(1)==ind(3)
        C = P_L(ind(1))*E_L(ind(2));
    end
end
function [E_W, V_W, cov_W] = original_scale(E_L,P_L)
% L_tr : 21 x 1 vector of cholesky space mean values
L_tr    = E_L;
L_tr_ii = L_tr(1:6);
P_tr  = diag(L_tr_ii);
P_tr(1,2) = L_tr(7); P_tr(1,3) = L_tr(8); P_tr(1,4) = L_tr(9); P_tr(1,5) = L_tr(10); P_tr(1,6) = L_tr(11);
P_tr(2,3) = L_tr(12); P_tr(2,4) = L_tr(13); P_tr(2,5) = L_tr(14); P_tr(2,6) = L_tr(15);
P_tr(3,4) = L_tr(16);P_tr(3,5) = L_tr(17); P_tr(3,6) = L_tr(18);
P_tr(4,5) = L_tr(19);P_tr(4,6) = L_tr(20);
P_tr(5,6) = L_tr(21);

cov_W  = P_tr*P_tr';
E_W_ii = diag(cov_W);
E_W_ij = [cov_W(1,2) cov_W(1,3) cov_W(1,4) cov_W(1,5) cov_W(1,6) cov_W(2,3) cov_W(2,4) cov_W(2,5) cov_W(2,6) cov_W(3,4) cov_W(3,5) cov_W(3,6)...
          cov_W(4,5) cov_W(4,6) cov_W(5,6)]';
E_W   =  [E_W_ii; E_W_ij];     
% Variance of the Covariance
V_tr    = diag(P_L);
V_tr    = diag(V_tr);
V_W     = diag(V_tr*V_tr');
end