addpath("quality_assess\","data\")
X = double(imread("mandril_color.tif")); % imshow(uint8(X));
X = double(imresize(X,[256,256])); 
R_1 = [10 10 3];
R_2 = R_1;

%% dtucker
maxiters = 50;
tol = 1e-15;
[~, ~, ~, ~, ~, Y] = dtucker_als(X, R_1, R_2,'maxiters',maxiters,'tol',tol);
psnr_ssim = quality_access(Y,X);
imshow(uint8(Y));
text(32,272,"28.xxxx/0.8xxx",'FontSize', 17)

num = numel(a)+numel(b);
for i = 1:size(cores,1)
    num = num + numel(cores{i});
end



%% Tucker hooi
addpath("tensor_toolbox-v3.2\")
X = tensor(X);
R = [80 80 3]; 
[T] = tucker_als(X,R);
G = T.core; U = T.u;
Y = full(T);
%e_hooi = norm(Y(:)-X(:))/norm(X(:));

num_tucker = numel(double(G));
for i = 1:size(U,1)
    num_tucker = num_tucker + numel(U{i});
end
MR4 = num_tucker/num;

%% tt and tw