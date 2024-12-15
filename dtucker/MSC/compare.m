sz = [11 15 15 11 3];
%% dtucker
r2 = [2 2 2 2 2]; r = r2;
num_dt = r(1)*sz(1)*r2(1) + r(2)*sz(2)*r2(2) + r(3)*sz(3)*r2(3) +...
    r(4)*sz(4)*r2(4) + r(5)*sz(5)*r2(5) + prod(r2)*2;

%% tucker
r = [2 9 4 7 2];
num_tucker = r(1)*sz(1) + r(2)*sz(2) + r(3)*sz(3) + r(4)*sz(4) +...
    + r(5)*sz(5) + prod(r);
re_tucker = num_tucker/num_dt;

%% TT
r = [2 7 9 3];
num_tt = sz(1)*r(1) + r(1)*sz(2)*r(2) + r(2)*sz(3)*r(3) +...
    r(3)*sz(4)*r(4) + r(4)*sz(5);
re_tt = num_tt/num_dt;

%% tw
r = [2 2 2 2 2]; l = [2 2 2 2 2];
num_tw = r(1)*sz(1)*l(1)*r(2) + r(2)*sz(2)*l(2)*(3) + r(3)*sz(3)*l(3)*r(4) +...
    r(4)*sz(4)*l(4)*(5) + r(5)*sz(5)*l(5)*(1) + prod(l);
re_tw = num_tw/num_dt;