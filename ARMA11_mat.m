function [M00,M10,M01,M02] = ARMA11_mat(y)

% Returns matrices M00, M10, M01, M11 for the rectangular MEP
%   (M00 + alfa1*M10 + gamma1*M01 + gamma1^2*M02) 
% whose eigenvalues are stationary points of the objective function for the 
% ARMA(1,1) model
 
class_t = superiorfloat(y);

N = length(y);
R = diag(ones(N-2,1,class_t),1) + diag(ones(N-2,1,class_t),-1);
ZB = zeros(N-1,class_t);
Id = eye(N-1,class_t);
zvec = zeros(N-1,1,class_t);
zrow = zeros(1,N-1,class_t);
y1 = y(1:N-1);
y2 = y(2:N);

M00 = [
    y2    Id    ZB    ZB;
    y1    ZB    Id    ZB;
    zvec  R     ZB    Id; 
    0     y1'   y2'   zrow;
    0     zrow  zrow  y2'
    ];

M10 = [
    y1    ZB    ZB    ZB;
    zvec  ZB    ZB    ZB;
    zvec  ZB    ZB    ZB; 
    0     zrow  y1'   zrow;
    0     zrow  zrow  y1'
    ];

M01 = [ 
    zvec  R     ZB    ZB;
    zvec  ZB    R     ZB;
    zvec  2*Id  ZB    R;
    0     zrow  zrow  zrow;
    0     zrow  zrow  zrow
    ];

M02 = [ 
    zvec  Id    ZB    ZB;
    zvec  ZB    Id    ZB;
    zvec  ZB    ZB    Id;
    0     zrow  zrow  zrow;
    0     zrow  zrow  zrow
    ];
