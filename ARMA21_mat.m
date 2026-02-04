function [M000,M100,M010,M001,M002] = ARMA21_mat(y)

% Returns matrices M000, M100, M010, M001, M002 for the rectangular MEP
%   (M000 + alfa1*M100 + alfa2*M010 + gamma1*M001 + gamma1^2*M002) 
% whose eigenvalues are stationary points of the objective function for the 
% ARMA(2,1) model

class_t = superiorfloat(y);

N = length(y);
R = diag(ones(N-3,1,class_t),1) + diag(ones(N-3,1,class_t),-1);
ZB = zeros(N-2,class_t);
Id = eye(N-2,class_t);
zvec = zeros(N-2,1,class_t);
zrow = zeros(1,N-2,class_t);
y1 = y(1:N-2);
y2 = y(2:N-1);
y3 = y(3:N);

M000 = [
    y3    Id    ZB    ZB    ZB;
    y2    ZB    Id    ZB    ZB;
    y1    ZB    ZB    Id    ZB;
    zvec  R     ZB    ZB    Id; 
    0     y2'   y3'   zrow  zrow;
    0     y1'   zrow  y3'   zrow;
    0     zrow  zrow  zrow  y3'
    ];

M100 = [
    y2    ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB; 
    zvec  ZB    ZB    ZB   ZB; 
    0     zrow  y2'   zrow zrow;
    0     zrow  zrow  y2'  zrow;
    0     zrow  zrow  zrow y2'
    ];

M010 = [
    y1    ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB; 
    zvec  ZB    ZB    ZB   ZB; 
    0     zrow  y1'   zrow zrow;
    0     zrow  zrow  y1'  zrow;
    0     zrow  zrow  zrow y1'
    ];

M001 = [ 
    zvec  R     ZB    ZB   ZB;
    zvec  ZB    R     ZB   ZB;
    zvec  ZB    ZB    R    ZB;
    zvec  2*Id  ZB    ZB    R;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow
    ];

M002 = [ 
    zvec  Id    ZB    ZB  ZB ;
    zvec  ZB    Id    ZB  ZB ;
    zvec  ZB    ZB    Id  ZB ;
    zvec  ZB    ZB    ZB  Id ;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow
    ];
