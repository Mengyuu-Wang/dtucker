function [C] = new_product(A,B)
    N = length(size(A));
    C = tensorprod(A,B,3,3);
    C = permute(C,[1 N 2 N+1 3:N-1]);
    if N >= 4
        szz = horzcat(size(A,1)*size(B,1),size(A,2)*size(B,2),size(A,4:N));
    else
        szz = horzcat(size(A,1)*size(B,1),size(A,2)*size(B,2));
    end
    C = reshape(C,szz);

end

%% A = randn([2 3 4 5 6]); B = randn([3 4 4]);