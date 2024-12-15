function[cores,a,b] = init_ele(sz,R_1,R_2)
%% [cores,a,b] = init_ele(sz,R_1,R_2);
    N = length(sz);
    a = randn(R_1);
    b = randn(R_2);
    cores = cell(N,1);
    for i = 1:N
        cores{i} = randn([R_1(i),sz(i),R_2(i)]);
    end
end