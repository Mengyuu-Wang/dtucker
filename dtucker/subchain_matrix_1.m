function G = subchain_matrix_1(a,b,cores,n)
    N = length(cores);
    if n == 1
        G = a;
        for i = 2:N
            G = tensorprod(G,cores{i},2,1);
        end
        G = tensorprod(G,b,2*(2:N)-1,2:N);
        G = reshape(G,size(a,1),numel(G)/(size(a,1)*size(b,1)),size(b,1));
    else
        G = a;
        for i = 1:n-1
            G = tensorprod(G,cores{i},1,1);
        end
        for i = n+1:N
            G = tensorprod(G,cores{i},2,1);
        end
        G = tensorprod(G,b,2*(2:N)-1,[1:n-1 n+1:N]);
        G = reshape(G,size(a,n),numel(G)/(size(a,n)*size(b,n)),size(b,n));
    end
    G = classical_mode_unfolding(G,2);
end