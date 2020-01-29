function rmat = hami2rmat(q)
%QUAT2LMAT Summary of this function goes here
%   Detailed explanation goes here



rmat = [ q(4),  q(3), -q(2), q(1);
        -q(3),  q(4),  q(1), q(2);
         q(2), -q(1),  q(4), q(3);
        -q(1), -q(2), -q(3), q(4)];

end