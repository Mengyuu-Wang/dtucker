function [a, b, cores, total_time, final_error, varargout] = dtucker_preSRFT_als(X, R_1, R_2, embedding_dims_1, embedding_dims_2, varargin)

iniTotalTime = tic;
%% Add relevant paths

addpath("mtimesx\mtimesx_20110223\");
addpath("tensor_toolbox-v3.2\");
mtimesx('SPEED');

% size of tensor
sz = size(X);

N = length(sz);
%% Handle inputs

% Optional inputs
params = inputParser;
addParameter(params, 'conv_crit', 'relative error');
addParameter(params, 'tol', 1e-3, @isscalar);
addParameter(params, 'maxiters', 50, @(x) isscalar(x) & x > 0);
addParameter(params, 'verbose', true, @isscalar);
parse(params, varargin{:});

conv_crit = params.Results.conv_crit;
tol = params.Results.tol;
maxiters = params.Results.maxiters;
verbose = params.Results.verbose;

[cores,a,b] = init_ele(sz,R_1,R_2);
re_cores = cell(N,1);
rre_cores = cell(N,1);
core_samples = cell(N,1);
core_samples_1 = cell(N,1);
core_samples_2 = cell(N,1);
for i = 1:N
    re_cores{i} = permute(cores{i},[2 3 1]);
end

%% Mixing tensor
diag_flips = cell(N,1);
for n = 1:N
    fft = dftmtx(sz(n))/sqrt(sz(n));
    diag_flips{n} = fft * diag((rand(sz(n),1)<0.5)*2-1);   
end

X_mixed = tensor(X);
for n = 1:N
    X_mixed = ttm(X_mixed, diag_flips{n},n);
end



sampling_probs = cell(N,1);
for n = 1:N
    sampling_probs{n} = ones(sz(n), 1)/sz(n);
end



% embedding_dims_1 = [100 100 100];
% embedding_dims_2 = [100 100 100];

breakup = ones(N,1);
slow_idx = cell(1,N);
sz_shifted = [1 sz(1:end-1)];
idx_prod = cumprod(sz_shifted);
sz_pts = cell(1,N);
for n = 1:N
    sz_pts{n} = round(linspace(0, sz(n), breakup(n)+1));
    slow_idx{n} = cell(1,breakup(n));
    for brk = 1:breakup(n)
        J = embedding_dims_1(n);
        samples_lin_idx_2 = prod(sz_shifted(1:n))*(sz_pts{n}(brk):sz_pts{n}(brk+1)-1).';
        slow_idx{n}{brk} = repelem(samples_lin_idx_2, J, 1);
    end
end

if nargout > 1 && tol > 0 && strcmp(conv_crit, 'relative error')
    conv_vec = zeros(1, maxiters);
    iter_time = zeros(1, maxiters+1);
end




%% Main loop
% Iterate until convergence, for a maximum of maxiters iterations

iter_time(1) = 0;
iniNopreTime = tic;
for it = 1:maxiters
    iniIterTime = tic;
    % Inner for loop
    for n = 1:N
        if n == 1
            % Construct sketch and sample cores
            J = embedding_dims_1(n);
            samples = nan(J, N);
            for m = 2:N
                samples(:,m) = randsample(sz(m), J, true, sampling_probs{m});
                core_samples{m} = re_cores{m}(samples(:,m), :, :);
            end
            a_re = permute(a,[2:N 1]);
            core_samples{1} = tensorprod(core_samples{2},a_re,3,1);
    
    
            % Compute the row rescaling factors
            rescaling = ones(J,1);
            for m = 2:N
                rescaling = rescaling ./ sqrt(sampling_probs{m}(samples(:,m)));
            end
            rescaling = rescaling ./ sqrt(J);
    
            % Construct sketched design matrix
            G_sketch = core_samples{1};
            for i = 3:N
                G_sketch = bridgeproduct(G_sketch,core_samples{i});
            end
            G_sketch = tensorprod(G_sketch,classical_mode_unfolding(b,1),2,2);
            G_sketch = classical_mode_unfolding(G_sketch,1);
            G_sketch = rescaling .* G_sketch;
    
    
            % Sample right hand side
            brk = 1;
            idx = [1:n-1 n+1:N];
            no_cols = sz_pts{n}(brk+1)-sz_pts{n}(brk);
            samples_lin_idx_1 = 1 + (samples(:, idx)-1) * idx_prod(idx).';
            samples_lin_idx = repmat(samples_lin_idx_1, no_cols, 1) + slow_idx{n}{brk};
            X_sampled = X_mixed(samples_lin_idx);
            Xn_sketch = reshape(X_sampled, J, no_cols);
    
    
            % Rescale right hand side
            Xn_sketch = rescaling .* Xn_sketch;

            Z = (G_sketch \ Xn_sketch).';

            cores{n} = getcore(Z, a, b, n); 
            re_cores{1} = permute(cores{1},[2 3 1]);
        else

            % Construct sketch and sample cores
            J = embedding_dims_1(n);
            samples = nan(J, N);
            for m = 1:N
                if m ~= n
                    samples(:, m) = randsample(sz(m), J, true, sampling_probs{m});
                    core_samples{m} = re_cores{m}(samples(:,m), :, :);
                end
            end
            a_re = permute(a,[1:n-1 n+1:N n]);
            core_samples{n} = tensorprod(core_samples{1},a_re,3,1);
    
            % Compute the row rescaling factors
            rescaling = ones(J, 1);
            for m = 1:N
                if m ~= n
                    rescaling = rescaling ./ sqrt(sampling_probs{m}(samples(:, m)));
                end
            end
            rescaling = rescaling ./ sqrt(J);
    
            % Construct sketched design matrix
            G_sketch = core_samples{n};
            for i = 2:N
                if i ~= n
                    G_sketch = bridgeproduct(G_sketch,core_samples{i});
                end
            end
            G_sketch = tensorprod(G_sketch,classical_mode_unfolding(b,n),2,2);
            G_sketch = classical_mode_unfolding(G_sketch,1);
            G_sketch = rescaling .* G_sketch;

    
    
            % Sample right hand side
            brk = 1;
            idx = [1:n-1 n+1:N];
            no_cols = sz_pts{n}(brk+1)-sz_pts{n}(brk);
            samples_lin_idx_1 = 1 + (samples(:, idx)-1) * idx_prod(idx).';
            samples_lin_idx = repmat(samples_lin_idx_1, no_cols, 1) + slow_idx{n}{brk};
            X_sampled = X_mixed(samples_lin_idx);
            Xn_sketch = reshape(X_sampled, J, no_cols);
    
    
            % Rescale right hand side
            Xn_sketch = rescaling .* Xn_sketch;

            Z = (G_sketch \ Xn_sketch).';
            cores{n} = getcore(Z, a, b, n); 
            re_cores{n} = permute(cores{n},[2 3 1]);
        end
    end
    %% updata A
    J_2 = embedding_dims_2(n);
    samples = nan(J_2, N);
    for i = 1:N
        rre_cores{i} = permute(cores{i},[2 1 3]);
    end
    for m = 1:N
        samples(:,m) = randsample(sz(m), J_2, true, sampling_probs{m});
        core_samples_1{m} = rre_cores{m}(samples(:,m), :, :);
    end
    
    % Compute the row rescaling factors
    rescaling = ones(J_2, 1);
    for m = 1:N
        rescaling = rescaling ./ sqrt(sampling_probs{m}(samples(:, m)));
    end
    rescaling = rescaling ./ sqrt(J_2);
    
    %%% 系数矩阵
    N_sketch = tensorprod(core_samples_1{1},b,3,1);
    for i = 2:N
        N_sketch = bridgeproduct(N_sketch,core_samples_1{i});
    end 
    N_sketch = rescaling .* N_sketch;

    % Sample right hand side
    X = tensor(X);
    Xn_sketch = X_mixed(samples);
    Xn_sketch = rescaling .* Xn_sketch;

        
    ZZ = (N_sketch \ Xn_sketch).';
    a = reshape(ZZ, size(a));

    %% updata B
    J_2 = embedding_dims_2(n);
    samples = nan(J_2, N);
    for m = 1:N
        samples(:,m) = randsample(sz(m), J_2, true, sampling_probs{m});
        core_samples_2{m} = re_cores{m}(samples(:,m), :, :);
    end
    
    % Compute the row rescaling factors
    rescaling = ones(J_2, 1);
    for m = 1:N
        rescaling = rescaling ./ sqrt(sampling_probs{m}(samples(:, m)));
    end
    rescaling = rescaling ./ sqrt(J_2);

    N_sketch = tensorprod(core_samples_2{1},a,3,1);
    for i = 2:N
        N_sketch = bridgeproduct(N_sketch,core_samples_2{i});
    end 
    N_sketch = rescaling .* N_sketch;

    % Sample right hand side
    X = tensor(X);
    Xn_sketch = X_mixed(samples);
    Xn_sketch = rescaling .* Xn_sketch;


    ZZZ = (N_sketch \ Xn_sketch).';
    b = reshape(ZZZ, size(b));


    iter_time(it+1) = iter_time(it) + toc(iniIterTime);
    for n = 1:N
        cores{n} = double(ttm(tensor(cores{n}),diag_flips{n}',2));
    end


    % Stop criterion: the original result also illustrates that one should avoid stopping too early when using ALS.
    %
    % Check convergence: Relative error
    if tol > 0 && strcmp(conv_crit, 'relative error')

        % Compute full tensor corresponding to cores
        Y = cores_2_tensor(a,b,cores,sz);

        % Compute current relative error
        er = norm(Y(:)-X(:))/norm(X(:));
        if verbose
            fprintf('\tRelative error after iteration %d: %.8f\n', it, er);
        end

        % Save current error to conv_vec if required
        if nargout > 1
            conv_vec(it) = er;
        end

        % Break if relative error below threshold
        if abs(er) < tol
            if verbose
                fprintf('\tRelative error below tol; terminating...\n');
            end
            break
        end


        % Check convergence: Norm change
        % Compute norm of TR tensor using normTR()
        % Code accompanying the paper "On algorithms for and computing with the tensor ring decomposition"
        % We delete this part because we will not use it as the stop criterion.
        % elseif


        % Just print iteration count
    else
        if verbose
            fprintf('\tIteration %d complete\n', it);
        end
    end

end
nopre_time = toc(iniNopreTime);
total_time = toc(iniTotalTime);
% final_error = er;
Y = cores_2_tensor(a,b,cores,sz);
final_error = norm(Y(:)-X(:))/norm(X(:));


if nargout > 3 && exist('conv_vec','var') && exist('iter_time','var')
    varargout{1} = conv_vec(1:it);
    varargout{2} = iter_time(2:it+1);
else
    varargout{1} = nan;
    varargout{2} = nan;
end

end
