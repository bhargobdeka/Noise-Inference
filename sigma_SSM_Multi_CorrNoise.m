clear;clc;
T=1000;                 %Time-serie length
n_x     = 6;
n_w     = n_x;
n       = n_x*2;
n_w2hat = n_x*(n_x+1)/2;
x_true  = zeros(n_x,T);      %Initialization of the vector of true values
y       = zeros(T,1);           %Initialization of the vector of observations
sV      = 0.001;
R       = sV^2;

%rng('default')
%% A Matrix
A_LL  = 1;
Block_A = cell(1,n_x);
for i = 1:n_x
    Block_A{i} = A_LL;
end
A = blkdiag(Block_A{:});
%% Q matrix
sW_LL = ones(1,n_x);
Block_Q = cell(1,n_x);
for i = 1:n_x
    Block_Q{i} = 2*sW_LL(i)^2;
end
% sW_LL1 = 1; sW_LL2 = 1; 
% sW_LL3 = 1; sW_LL4 = 1;
% Q_LL    = {sW_LL1^2 sW_LL2^2 sW_LL3^2 sW_LL4^2};
Q       = blkdiag(Block_Q{:});
true_value = 0.5;
% % % Off-Diagonal terms in Q
% sW_AR12 = true_value;
% sW_AR13 = true_value;
% sW_AR23 = true_value;
% sW_AR14 = true_value; sW_AR24 = true_value; sW_AR34 = true_value ;
% 
% Q(1,2)  = sW_AR12; Q(2,1)  = sW_AR12;
% Q(1,3) = sW_AR13; Q(2,3) = sW_AR23; Q(3,1) = sW_AR13; Q(3,2) = sW_AR23;
% Q(1,2:4) = [sW_AR12 sW_AR13 sW_AR14];
% Q(2:4,1) = [sW_AR12 sW_AR13 sW_AR14]';
% Q(2,3:4) = [sW_AR23 sW_AR24];
% Q(3:4,2) = [sW_AR23 sW_AR24];
% Q(3,4)   = sW_AR34;
% Q(4,3)   = sW_AR34;

% % Off-Diagonal terms in Q
sW_AR12 = true_value;
sW_AR13 = true_value;
sW_AR14 = true_value;
sW_AR15 = true_value;
sW_AR16 = true_value;
sW_AR23 = true_value;
sW_AR24 = true_value; 
sW_AR25 = true_value;
sW_AR26 = true_value; 
sW_AR34 = true_value;
sW_AR35 = true_value;
sW_AR36 = true_value;
sW_AR45 = true_value;
sW_AR46 = true_value;
sW_AR56 = true_value;

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
EX=zeros(n_x+n_w+n_w2hat,T);
E_w2hat = [[5.3100e+00   5.1600e+00   5.0100e+00   4.8600e+00   4.7100e+00   4.5600e+00 ] zeros(1,n_w2hat-n_x)];
%P_w2hat = 0.5*ones(1,n_w2hat);
P_w2hat = [2.8260e+00   2.8710e+00   2.8260e+00   2.7810e+00   2.7360e+00   2.6910e+00 ...
           0.774 7.5150e-01   7.2900e-01   7.0650e-01   6.8400e-01   7.5150e-01   7.2900e-01   7.0650e-01 ...
           6.8400e-01   7.2900e-01   7.0650e-01   6.8400e-01   7.0650e-01   6.8400e-01   6.8400e-01];
EX(:,1)=[zeros(1,n_x) nan*zeros(1,n_x) E_w2hat]';                      
PX(:,:,1)=diag([0.1*ones(1,n_x),nan*zeros(1,n_x),P_w2hat]);               
%my   = zeros(T,2);
% SY  = zeros(T,2);

%% Prediction 
for t=2:T
    Ep           = [A*EX(1:n_x,t-1);zeros(n_x,1)];            % mu_t|t-1
    s_w_sq       = EX(end-n_w2hat+1:end,t-1);
    n_V          = n_x;    n_Cov = n_w2hat-n_x;               % no of variance terms and covariance terms for W2hat
    s_wii        = s_w_sq(1:n_V);                             % variance terms 
    s_wij        = s_w_sq(n_V+1:end);                         % covariance terms
    Sx           = A*PX(1:n_x,1:n_x,t-1)*A';         
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
%Check PSD
% try
%     chol(P_W);
% catch
%     warning('Not a PSD matrix');
%     P_W = nearestSPD(P_W);
% end

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
        % Check PSD
%         try 
%             chol(PX_wy);
%         catch
%             warning('Not a PSD matrix');
%             PX_wy  = (PX_wy+PX_wy')/2;
%             PX_wy  = nearestSPD(PX_wy);
%         end
        
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
        %Check PSD
%         try 
%             chol(PX_wpy);
%         catch
%             warning('Not a PSD matrix');
%             PX_wpy               = (PX_wpy + PX_wpy')/2;
%             PX_wpy               = nearestSPD(PX_wpy);
%         end
%         try 
%             chol(PX_wpy);
%         catch
%             warning('Not a PSD matrix');
%         end
        
        %% Creating E[W^p]
        EX_wpy = [m_wii_y' m_wiwj_y];

%% Computing prior mean and covariance matrix of Wp
        m_wsqhat    = EX(end-n_w2+1:end,t-1);
        s_wsqhat    = PX(end-n_w2+1:end,end-n_w2+1:end,t-1);
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
        cov_wp             = cell2mat(reshape(cov_prior_wijkl,size(cov_prior_wijkl,2),1));
        s_wp                 = zeros(size(PX_wp,1));
        s_wp(1:end-1,2:end)  = cov_wp;
        PX_wp               = PX_wp + s_wp;
        PX_wp               = triu(PX_wp)+triu(PX_wp,1)'; % adding the lower triangular matrix
%         try 
%             chol(PX_wp);
%         catch
%             warning('Not a PSD matrix');
%             PX_wp                = (PX_wp+PX_wp')/2;
%             PX_wp                = nearestSPD(PX_wp);
%         end
%         %Check PSD
%         try 
%             chol(PX_wp);
%         catch
%             warning('Not a PSD matrix');
%         end
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
        EX(end-n_w2+1:end,t)                 = EX(end-n_w2+1:end,t-1)                  +      J*(E_Wp_pos' - E_Wp_prior);
        PX(end-n_w2+1:end,end-n_w2+1:end,t)  = PX(end-n_w2+1:end,end-n_w2+1:end,t-1)   +      J*(P_Wp_pos - P_Wp_prior)*J';
        PX(:,:,t)                            = nearestSPD(PX(:,:,t));
        
    end
end
varii = EX(end-n_w2+1:end-n_w2+n_x,T);
varij = EX(end-n_w2+n_x+1:end,T);
estim_sigma_w = [sqrt(varii);varij];

%% Plotting
%  Plotting Variances
t  = 1:length(EX);
figure;
for i=1:n_x
    subplot(3,2,i)
    xw = EX(n+i,t);
    sX = permute(sqrt(PX(n+i,n+i,t)),[1,3,2]);
    plot(t,repmat(2*sW_LL(i)^2,[1,length(EX)]),'-.r','Linewidth',1)
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
    xw = EX(n+n_w+i,t);
    sX = permute(sqrt(PX(n+n_w+i,n+n_w+i,t)),[1,3,2]);
    plot(t,repmat(true_value,[1,length(EX)]),'-.r','Linewidth',1)
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