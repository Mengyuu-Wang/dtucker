function [S,history]=MSC(X,paras)

lambda=paras.lambda;
rX=paras.rX;
R_mera=paras.R;
% Initialize...

K = length(X); N = size(X{1},2); %sample number


for k=1:K
    Z{k} = zeros(N,N); %Z{2} = zeros(N,N);
    Gamma{k} = zeros(N,N);
    Y{k} = zeros(N,N);
    E{k} = zeros(size(X{k},1),N); %E{2} = zeros(size(X{k},1),N);
    Lambda{k} = zeros(size(X{k},1),N); %Y{2} = zeros(size(X{k},1),N);
        L(k)=size(X{k},1);
end

w = zeros(N*N*K,1);
g = zeros(N*N*K,1);
dim1 = N;dim2 = N;dim3 = K;
sX = [N, N, K];

%set Default
parOP         =    false;
ABSTOL        =    1e-6;
RELTOL        =    1e-4;

Isconverg = 0;epson = 1e-6;

iter = 0;
mu_2 = 5e-4; max_mu_2 = 10e10; pho_mu_2 = 2;
mu_1 = 10e-4; max_mu_1 = 10e10; pho_mu_1 = 2;
t1 = cputime;
% tic;

while(Isconverg == 0)
      iter = iter + 1;
    fprintf('processing iter %d\n', iter);
    for k=1:K
        %1 update Z^k
        tmp = X{k}'*Lambda{k} + mu_2*X{k}'*X{k} - mu_2*X{k}'*E{k}  +  mu_1*Y{k}- Gamma{k};
        Z{k}=inv(mu_1*eye(N,N)+ (mu_2)*X{k}'*X{k})*tmp;
        E{k}=solve_l1l2(X{k}-X{k}*Z{k}+Lambda{k}/mu_2,lambda/mu_2);
        Lambda{k} = Lambda{k} + mu_2*(X{k}-X{k}*Z{k}-E{k});

    end
    Z_tensor = cat(3, Z{:,:});
    G_tensor = cat(3, Gamma{:,:});
    z = Z_tensor(:);
    g = G_tensor(:);
    Y_tensor=Z_tensor;
    
    ZZ = reshape(Z_tensor + 1/mu_1*G_tensor,rX);
    
    maxiters = 50;
    tol = 1e-15;

    %% dtucker
    [~, ~, ~, ~, ~, YYY] = dtucker_als(ZZ,R_mera{1},R_mera{1});

    %% Tukcer
    % [T] = tucker_als(tensor(ZZ),R_mera{1});
    % YYY = double(full(T));

    %% TT
    % X_hat_arr = TTOI(ZZ,R_mera{1},maxiters,tol);
    % YYY = X_hat_arr{3};

    %% TW
    % [~,~,~,~,YYY] = TW_ALS(ZZ, R_mera{1}, [2 2 2 2 2],'tol',tol, 'maxiters', maxiters);
     
    Y_tensor = reshape(YYY,sX);
    y = Y_tensor(:);
    g = g + mu_1*(z - y);


    %% coverge condition
%     Isconverg = 1;
    for k=1:K
            norm_Z(k) = norm(X{k}-X{k}*Z{k}-E{k},inf);
%             fprintf('  RE %7.10f   \n  ', norm_Z(k) );

     
        Y{k} = Y_tensor(:,:,k);
        G_tensor = reshape(g, sX);
        Gamma{k} = G_tensor(:,:,k);
        
  
            norm_Z_G(k)  = norm(Z{k}-Y{k},inf);
%             fprintf('ME %7.10f    \n', norm_Z_G(k) );
    
    end
    

      history.re(iter)=max(norm_Z);
      history.me(iter)=max(norm_Z_G);
     
       fprintf('RE %10.16f    \n', history.re(iter));
       
    if history.re(iter) <epson && history.me(iter)< epson
         Isconverg  = 1;
    end
    
    if (iter>50)
        Isconverg  = 1;
    end
  
    mu_2 = min(mu_2*pho_mu_2, max_mu_2);
    mu_1 = min(mu_1*pho_mu_1, max_mu_1);
end


S = 0;
for k=1:K
    S = S + abs(Z{k})+abs(Z{k}');
end

