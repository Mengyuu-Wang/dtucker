function [C] = bridgeproduct(A,B)
    sz = size(A);
    N = length(sz);
    A = permute(A,[2 4:N 3 1]);
    A = reshape(A,numel(A)/(sz(1)*sz(3)),sz(3),sz(1));
    B = permute(B,[3 2 1]);
    C = mtimesx(A,B);
    szz = horzcat(sz(2),prod(sz(4:N)),size(B,2),size(B,3));
    C = permute(reshape(C,szz),[4,1,3,2]);
    szz = horzcat(size(B,3),sz(2)*size(B,2),sz(4:N));
    C = reshape(C,szz);
end

%% A = randn([6 5 7 8 9]); B = randn([6 3 7]);
