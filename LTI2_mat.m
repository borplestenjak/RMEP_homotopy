function [M00,M10,M01,M20,M11,M02] = LTI2_mat(y)

% [M00,M10,M01,M20,M11,M02] = LTI_matrices(y) returns (3N-4)x(3N-5) 
% matrices M00, M10, M01, M20, M11, M02 for the rectangular MEP
% (M00 + l*M10 + u*M01 + l^2*M20 + l*u*M11 + u^2*M02)z = 0
% for the optimal LTI(2) identification problem

class_t = superiorfloat(y);

N = length(y);
y = y(:);
y1 = y(1:N-2);
y2 = y(2:N-1);
y3 = y(3:N);
Id = eye(N-2,class_t);
ZM = zeros(N-2,class_t);
Zr = zeros(1,N-2,class_t);
Zc = zeros(N-2,1,class_t);

P010 = diag(ones(N-3,1,class_t),1) + diag(ones(N-3,1,class_t),-1);
if N==3
    P100 = numeric_t(0,class_t);
else
    P100 = diag(ones(N-4,1,class_t),2) + diag(ones(N-4,1,class_t),-2);
end

M00 = [y3   Id     ZM    ZM
       y2   P010   Id    ZM
       y1   P100   ZM    Id
       0    y2.'   y3.'  Zr
       0    y1.'   Zr    y3.'];

M10 = [y2   P010   ZM    ZM
       Zc   2*Id   P010  ZM
       Zc   P010   ZM    P010
       0    Zr     y2.'  Zr
       0    Zr     Zr    y2.'];

M01 = [y1   P100   ZM    ZM
       Zc   P010   P100  ZM
       Zc   2*Id   ZM    P100
       0    Zr     y1.'  Zr
       0    Zr     Zr    y1.'];

M11 = [Zc   P010   ZM    ZM
       Zc   ZM     P010  ZM
       Zc   ZM     ZM    P010
       0    Zr     Zr    Zr
       0    Zr     Zr    Zr];

M20 = [Zc   Id     ZM    ZM
       Zc   ZM     Id    ZM
       Zc   ZM     ZM    Id
       0    Zr     Zr    Zr
       0    Zr     Zr    Zr];

M02 = M20;

end

function e = LTI2_err(y,alpha1,alpha2,class_t)

    N = length(y);

    tmp1 = [alpha2*eye(N-2,class_t) zeros(N-2,2,class_t)];
    tmp2 = [zeros(N-2,1,class_t) alpha1*eye(N-2,class_t) zeros(N-2,1,class_t)];
    tmp3 = [zeros(N-2,2,class_t) eye(N-2,class_t)];
    TA = tmp1 + tmp2 + tmp3;
    DC = TA*TA';

    err = TA'*(DC\(TA*y));
    e = norm(err)^2;

end