function varargout = RGFFCM(X, K, V, G, varargin)
% RGFFCM - Fuzzy C-Means Clustering Based on Regional Feature Grouping.
%
% Syntax
% =================
% out = RGFFCM(X, K, V, G, 'c_idx', [], 'debug_mode', 1)
%
% Input Arguments
% =================
% X             Data to be clustered, NxM matrix.
% K             Number of clusters.
% V             The weights of the feature groups.
% G             The labels of feature groups. 
%
% c_idx         The index of initial cluster centers, the default is [].
% version       Implemented code version, {'CGFFCM', } 
% debug_mode    Control debug information, 0: Silent, 1: Call information, 2: Call details.	
%
% Output Arguments
% =================
% out           L or {U, C} or {U, C, W, Z}.
%
% References
% =================
% [1] A.Golzari oskouei, M.Hashemzadeh, B.Asheghi  and M.Balafar, 
% "CGFFCM: Cluster-weight and Group-local Feature-weight learning in 
% Fuzzy C-Means clustering algorithm for color image segmentation", 
% Applied Soft Computing, vol. 113, pp. 108005, 2021.
% [2] M. Hashemzadeh, A. Golzari Oskouei, and N. Farajzadeh, 
% "New fuzzy C-means clustering method based on feature-weight and cluster-weight learning," 
% Appl. Soft Comput., vol. 78, pp. 324-345, 2019/05/01/, 2019.
% [3] https://github.com/Amin-Golzari-Oskouei/CGFFCM
% [4] https://github.com/Amin-Golzari-Oskouei/FWCW_FCM


% Parameter Initialization    
% =========================================================
p = inputParser; fun_name = 'RGFFCM';  
% addParameter(p,'V', [0.2 0.6 0.2]); 
% addParameter(p,'G', [1 1 1 2 2 2 3 3 3]); 
addParameter(p,'c_idx',[]);
addParameter(p,'t_max',100);
addParameter(p,'p_init',0);
addParameter(p,'p_max',0.5); 
addParameter(p,'p_step',0.01); 
addParameter(p,'eps_',1e-6); 
addParameter(p,'q',-4); 
addParameter(p,'alpha_',2); 
addParameter(p,'beta_',0.3); 
addParameter(p,'I', 0.0001); 
addParameter(p,'version','CGFFCM');
addParameter(p,'debug_mode',0); 
parse(p,varargin{:});

if p.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name);
elseif p.Results.debug_mode == 2     
    fprintf('\nCall functions:\t%s\n', fun_name);
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(p.Results);
end

if isempty(p.Results.c_idx)
    rng(123); c_idx = randperm(N, K);
else
    c_idx = p.Results.c_idx;
end


% Method Implementation 
% =========================================================
if strcmpi(p.Results.version, 'CGFFCM')
    [U,C,W,Z,F_history,U_history,W_history,Z_history] = RGFFCM_v1(X,K,V,G,c_idx,...
        p.Results.t_max,p.Results.p_init,p.Results.p_max,p.Results.p_step,...
        p.Results.eps_,p.Results.q,p.Results.alpha_,p.Results.beta_,...
        p.Results.I);  % ,p.Results.V,p.Results.G
end

% Output Settings
% =========================================================
if nargout == 8
    varargout = {U,C,W,Z,F_history,U_history,W_history,Z_history};
elseif nargout == 5
    varargout = {U,C,W,Z,F_history};
elseif nargout == 4
    varargout = {U,C,W,Z};
elseif nargout == 2
    varargout = {U,C};
else
    [~, idx] = max(U,[],2);
    varargout = {idx};
end   

end

%% Core Function
function [U,C,W,Z,F_history,U_history,W_history,Z_history] = RGFFCM_v1(X,K,V,G,c_idx,t_max,p_init,p_max,p_step,eps_,q,alpha_,beta_,I)
%
% This demo implements the CGFFCM algorithm as described in
% [1] A.Golzari oskouei, M.Hashemzadeh, B.Asheghi  and M.Balafar, "CGFFCM: Cluster-weight 
% and Group-local Feature-weight learning in Fuzzy C-Means clustering algorithm for color
% image segmentation", Applied Soft Computing, vol. 113, pp. 108005, 2021.
% [2] M. Hashemzadeh, A. Golzari Oskouei, and N. Farajzadeh, 
% "New fuzzy C-means clustering method based on feature-weight and cluster-weight learning," 
% Appl. Soft Comput., vol. 78, pp. 324-345, 2019/05/01/, 2019.
% 
% Arguments:
% ===============
% X         an NxM data matrix, where each row corresponds to an instance.
% K         the number of clusters.
% V         the weights of feature groups (G).
% G         the labels of feature groups. 
% c_idx     Index of initial cluster centers
% t_max 	the maximum number of iterations. 
% p_init 	the initial p value (0<=p_init<1).
% p_max 	the maximum admissible p value (0<=p_max<1).
% p_step 	the step for increasing p (p_step>=0).
% eps_      the tolerance of the objective function.
% q         the power value of the feature weight(in paper).
% alpha_    the fuzzy degree (fuzzification coefficient).
% beta_     amount of memory for the weights updates.
% I         used to calculate gamma_m, value interval (0,1].

% 
% Returns
% ================
% U             the membership matrix (N,K)
% C             the cluster centers matrix(K,M)
% W             the feature weight matrix (K,M)
% Z             the cluster weight vector (1,K)
% F_history     records the objective function (1,t)
% U_history     records the membership matrix (N,tK)
% W_history     records the feature weight (tK,M)
% Z_history     records the cluster weight (t,K)
%--------------------------------------------------------------------------

if p_init<0 || p_init>=1
    error('p_init must take a value in [0,1)');
end

if p_max<0 || p_max>=1
    error('p_max must take a value in [0,1)');
end

if p_max<p_init
    error('p_max must be greater or equal to p_init');
end

if p_step<0
    error('p_step must be a non-negative number');
end

if beta_<0 || beta_>1
    error('beta must take a value in [0,1]');
end

if q==0
    error('beta must be a non-zero number');
end

if p_init==p_max    
    if p_step~=0
        fprintf('p_step reset to zero, since p_max equals p_init\n\n');
    end    
    p_flag=0;
    p_step=0;  
elseif p_step==0    
    if p_init~=p_max
        fprintf('p_max reset to equal p_init, since p_step=0\n\n');
    end    
    p_flag=0;
    p_max=p_init;    
else
    p_flag=1; % p_flag indicates whether p will be increased during the iterations.
end

if min(G)~=1 || max(G)~=length(unique(G))
    error('G must be from 1 to N');
end

%--------------------------------------------------------------------------

% Initialization parameters
t=0;                    % Number of iterations
% p_init=0;             % Initial p_init.
empty=0;                % Count the number of iterations for which an empty or singleton cluster is detected.
p=p_init;               % Initial p value.
p_prev=p-10^(-8);       % Dummy value.
% eps_=1e-6;            % Tolerance of the objective function.

% Initialize cluster centers
[N, M] = size(X); 
if isempty(c_idx)
    rng(123); c_idx = randperm(N, K);
end
C = X(c_idx, :);

% Initialize weights
Z = ones(1,K)/K; W = zeros(K,M); v = zeros(1,M);         
if (length(G)==M) && (length(unique(G))==length(V))    
    for i=1:length(unique(G))
        W(:,G==i) = 1/nnz(G==i); v(G==i)=V(i);
    end
else
    error('Weights of features do not match');
end

% Iterative calculation preparation
% I = 0.0001;           % interval (0,1]
gamma_m = I ./ var(X);  % the inverse variance of the m-th feature
gamma_m(isinf(gamma_m)) = max(gamma_m(~isinf(gamma_m)))+1;
 
% Other initializations.
F_old=inf;              % Previous iteration objective (used to check convergence).
U_history=[]; Z_history=[]; W_history=[]; F_history=[];

D_exp = zeros(N,M,K); D_nmk = zeros(N,M,K); D_nk = zeros(N,K); 
D_omega = zeros(K,M); D_z = zeros(1,K);

%--------------------------------------------------------------------------

% The CGFFCM iterative procedure.
fprintf('\nStart of RFGFCM iterations\n');
fprintf('----------------------------------\n\n');

while 1
    
    t=t+1; 
    
    % Calculate non-Euclidean distance.
    omega_q = transpose(W.^q); omega_q(isinf(omega_q)) = 0;     % (M,K)
    for k=1:K
        D_exp(:,:,k) = exp((-1.*repmat(gamma_m,N,1)).*((X-repmat(C(k,:),N,1)).^2));
        D_nmk(:,:,k) = (1-D_exp(:,:,k)).*repmat(v,N,1);
        D_nk(:,k) = Z(1,k).^p * D_nmk(:,:,k) * omega_q(:,k);
    end
    
    % Update the cluster assignments matrix U.  
    D_u = zeros(N,K);
    for k=1:K
        tmp = (D_nk./repmat(D_nk(:,k),1,K)).^(1/(alpha_-1));  %????
        tmp(isinf(tmp))=0; tmp(isnan(tmp))=0;
        D_u = D_u+tmp;
    end  
    
    U=transpose(1./D_u); U(isnan(U))=1; U(isinf(U))=1;   % (K,N)     
    if nnz(D_nk==0)>0
        for i=1:N
            if nnz(D_nk(i,:)==0)>0
                U(find(D_nk(i,:)==0),i) = 1/nnz(D_nk(i,:)==0);
                U(find(D_nk(i,:)~=0),i) = 0;
            end
        end
    end
    
    % Calculate the CGFFCM objective.
    U_alpha = U.^alpha_;   % (K,N) 
    F = sum(sum(D_nk.*transpose(U_alpha)));
    F_history(t) = F;
    
    % If empty or singleton clusters are detected after the update.
    for k=1:K
        tmp = find(U(k,:)<=0.05);
        if length(tmp)==N-1 || length(tmp)==N
            
            fprintf('Empty or singleton clusters detected for p=%g.\n',p);
            fprintf('Reverting to previous p value.\n\n');
            
            F = NaN; % Current objective undefined.
            empty=empty+1;
            
            % Reduce p when empty or singleton clusters are detected.
            if empty>1
                p=p-p_step;              
                % The last p increase may not correspond to a complete p_step,
                % if the difference p_max-p_init is not an exact multiple of p_step.
            else
                p=p_prev;
            end
            
            p_flag=0; %Never increase p again.
            
            % p is not allowed to drop out of the given range.
            if p<p_init || p_step==0
                
                fprintf('\n+++++++++++++++++++++++++++++++++++++++++\n\n');
                fprintf('p cannot be decreased further.\n');
                fprintf('Either p_step=0 or p_init already reached.\n');
                fprintf('Aborting Execution\n');
                fprintf('\n+++++++++++++++++++++++++++++++++++++++++\n\n');
                
                % Return NaN to indicate that no solution is produced.
                C = NaN(K,M); 
                U = transpose(U); U_history = transpose(U_history);
                return;
            end
            
            % Continue from the assignments and the weights corresponding to the decreased p value.
            a=(K*empty)-(K-1); b=K*empty;
            U=U_history(a:b,:);
            Z=Z_history(empty,:);           
            W=W_history(a:b,:);
            break;
        end
    end
    
    if ~isnan(F) && mod(t,5)==0
        fprintf('p=%g\n',p);
        fprintf('The RGFFCM objective is F=%f\n\n',F);
    end
    
    % Check for convergence. Never converge if in the current (or previous)
    % iteration empty or singleton clusters were detected.
    if ~isnan(F) && ~isnan(F_old) && (abs(1-F/F_old)<eps_ || t>=t_max)        
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
        fprintf('Converging for p=%g after %d iterations.\n',p,t);
        fprintf('The final RGFFCM objective is F=%f.\n',F);
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');       
        break;       
    end    
    F_old=F;
    
    % Update the cluster centers.   
    for k=1:K
        C(k,:) = (U_alpha(k,:)*(X.*D_exp(:,:,k)))./(U_alpha(k,:)*D_exp(:,:,k));
    end
    
    %Increase the p value.
    if p_flag==1        
        % Keep the assignments-weights corresponding to the current p.
        % These are needed when empty or singleton clusters are found in
        % subsequent iterations.
        U_history = [U; U_history];
        W_history = [W; W_history];
        Z_history = [Z; Z_history];

        p_prev=p; p=p+p_step;        
        if p>=p_max
            p=p_max; p_flag=0;
            fprintf('p_max reached\n\n');
        end
    end
    W_old=W; Z_old=Z;
    
    % Update the feature weight matrix W.
    for k=1:K
        D_omega(k,:) = U_alpha(k,:) * (1 - D_exp(:,:,k));
    end
    
    W = zeros(K,M); G_labels = unique(G);
    for i=1:length(G_labels)
        g = G_labels(i); n_g = nnz(G==g);
        tmp1 = zeros(K, n_g); D_w = D_omega(:,G==g);
        for j=1:n_g
            tmp2 = (D_w./repmat(D_w(:,j),1,n_g)).^(1/(q-1));
            tmp2(isnan(tmp2))=0; tmp2(isinf(tmp2))=0;
            tmp1 = tmp1+tmp2;
        end
        W(:,G==g) = 1./tmp1;
    end
    W(isnan(W))=1; W(isinf(W))=1;   
    
    if nnz(D_omega==0)>0
        for k=1:K
            for j=1:length(unique(G))
                g = G_labels(j);
                if nnz(D_omega(k,G==g)==0)>0
                    W(k,find(G==g & D_omega(k,:)==0)) = 1/nnz(D_omega(k,G==g)==0);
                    W(k,find(G==g & D_omega(k,:)~=0)) = 0;
                end
            end
        end
    end

    
    % Update the cluster weights.   
    for k=1:K
        D_z(1,k) = U_alpha(k,:) * D_nmk(:,:,k) * omega_q(:,k);           
    end 
    
    tmp = sum((repmat(D_z,K,1)./transpose(repmat(D_z,K,1))).^(1/(p-1)));
    tmp(isinf(tmp))=0; tmp(isnan(tmp))=0;
    
    Z = 1./tmp; Z(isnan(Z))=1; Z(isinf(Z))=1;    
    if nnz(D_z==0)>0
        Z(D_z==0) = 1/nnz(D_z==0);
        Z(D_z~=0) = 0;
    end    
    
    % Memory effect.
    W = (1-beta_)*W + beta_*W_old;
    Z = (1-beta_)*Z + beta_*Z_old;
end

U = transpose(U); U_history = transpose(U_history);

end