function[a, b, cores, total_time, final_error, Y] = dtucker_als(X, R_1, R_2, varargin)
    
    sz = size(X);

    params = inputParser;
    addParameter(params, 'conv_crit', 'relative error');
    addParameter(params, 'tol', 1e-15, @isscalar);
    addParameter(params, 'maxiters', 50, @(x) isscalar(x) & x > 0);
    addParameter(params, 'verbose', true, @isscalar);
    addParameter(params, 'random_init', true, @isscalar);
    parse(params, varargin{:});
    
    conv_crit = params.Results.conv_crit;
    tol = params.Results.tol;
    maxiters = params.Results.maxiters;
    verbose = params.Results.verbose;
    random_init = params.Results.random_init;
    
    if random_init == true
        [cores,a,b] = init_ele(sz,R_1,R_2);
    else
        [a, b, cores] = dtucker_svd(X, R_1, R_2);
    end
    
    iniTotalTime = tic;

    % nargout： the output variables of Function
    if nargout > 1 && tol > 0 && strcmp(conv_crit, 'relative error')
        conv_vec = zeros(1, maxiters);
        iter_time = zeros(1, maxiters+1);
    end

    N = length(sz);
    iter_time(1) = 0;
    for it = 1:maxiters
        iniIterTime = tic;
        for n = 1:N
            G = subchain_matrix_1(a, b, cores, n);
                        
            XnT = classical_mode_unfolding(X, n).';
            
            Z = (G \ XnT).';
            cores{n} = getcore(Z,a,b,n);        
        end
        %% updata A
        C = subchain_matrix_2(b,cores,sz).';
        
        XT = reshape(X,numel(X),1);
        
        ZZ = (C \ XT).';
        a = reshape(ZZ, size(a));

        %% updata B
        C = subchain_matrix_3(a,cores,sz);
        
        %XT = reshape(X,numel(X),1);
        
        ZZ = (C \ XT).';
        b = reshape(ZZ, size(b));

        
        iter_time(it+1) = iter_time(it) + toc(iniIterTime);
    
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






