function core = getcore(Z,a,b,n)
    core = permute(reshape(Z,numel(Z)/(size(a,n)*size(b,n)),size(a,n),size(b,n)),[2 1 3]);
end