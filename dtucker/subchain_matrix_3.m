function C = subchain_matrix_3(a,cores,sz)
    N = length(cores);
    C = a;
    for i = 1:N
        C = tensorprod(C,cores{i},1,1);
    end
    C = permute(C,[1:2:2*N 2:2:2*N]);
    C = reshape(C,prod(sz),numel(C)/prod(sz));
end


%% test
function result = tensor_transformation(cores, a)
    % tensors: a cell array of N 3D tensors, where each tensor has dimensions Js x Is x Ls
    % big_tensor: an N-dimensional tensor of size J_1 x J_2 x ... x J_N

    N = numel(cores);
    Is_dims = zeros(1, N);
    Ls_dims = zeros(1, N);
    Js_dims = zeros(1, N);

    % Extract the Is, Js, and Ls dimensions from each tensor
    for s = 1:N
        [Js, Is, Ls] = size(cores{s});
        Is_dims(s) = Is;
        Ls_dims(s) = Ls;
        Js_dims(s) = Js;
    end

    % Reshape big_tensor to ensure it matches the Js dimensions
    big_tensor_reshaped = reshape(a, Js_dims);

    % Initialize the result tensor with big_tensor_reshaped
    result = big_tensor_reshaped;
    result = tensor(result);
    % Perform the tensor contraction to remove Js and combine Is and Ls
    for s = 1:N
        tensor_s = cores{s};

        % Reshape the tensor_s to (Js_s) x (Is_s * Ls_s)
        tensor_s = reshape(tensor_s, Js_dims(s), Is_dims(s) * Ls_dims(s));

        % Tensor contraction along the Js dimension
        result = ttm(result, tensor_s', s); % Multiply along the s-th dimension
    end


% Permute dimensions to group Is and Ls together
perm_order = reshape(1:2*N, 2, []).';
perm_order = perm_order(:).';

result = permute(result, perm_order);
end
