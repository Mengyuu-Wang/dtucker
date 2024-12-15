function C = subchain_matrix_2(b,cores,sz)
    N = length(cores);
    C = b;
    for i = 1:N
        C = tensorprod(C,cores{i},1,3);
    end
    C = permute(C,[1:2:2*N 2:2:2*N]);
    C = reshape(C,numel(C)/prod(sz),prod(sz));
end