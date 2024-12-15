function Y = cores_2_tensor(a,b,cores,sz)
    C = subchain_matrix_2(b,cores,sz);
    a = reshape(a,numel(a),1);
    Y = tensorprod(C,a,1,1);
    Y = reshape(Y,sz);
end