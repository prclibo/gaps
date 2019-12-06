function d = detv(B)
% MULTIPOL determinant function. This is just a quick hack implementing
% determinants of size up to 2,3,4, 5, 6.

nn = size(B,1);

if nn == 2
    d = B(1) * B(4) - B(2) * B(3);
elseif nn == 3
    cp = [2 3;
        1 3 ;
        1 2 ;
        ];
    rr = [[1 2 3 ]'];
    
    d = multipol();
    % d = 0;
    for i = 1:3
        subd(i) = detv(B(2:end,cp(i,:)));
        d       = (-1)^(rr(i)-1) * B(1,rr(i))*subd(i) + d;
    end
elseif nn == 4;
    cp = [2 3 4;
        1 3 4;
        1 2 4;
        1 2 3];
    rr = [[1 2 3 4]'];
    
    d = multipol();
    % d = 0;
    for i = 1:4
        subd(i) = detv(B(2:end,cp(i,:)));
        d       = (-1)^(rr(i)-1) * B(1,rr(i))*subd(i) + d;
    end
elseif nn == 5
    cp = [ 2 3 4 5;
        1 3 4 5;
        1 2 4 5;
        1 2 3 5;
        1 2 3 4];
    rr = [[1 2 3 4 5]'];
    
    d = multipol();
    % d = 0;
    for i = 1:5
        subd(i) = detv(B(2:end,cp(i,:)));
        d       = (-1)^(rr(i)-1) * B(1,rr(i))*subd(i) + d;
    end
elseif nn == 6;
    cp = [ 2 3 4 5 6;
        1 3 4 5 6;
        1 2 4 5 6;
        1 2 3 5 6;
        1 2 3 4 6;
        1 2 3 4 5];
    rr = [[1 2 3 4 5 6]'];
    
    d = multipol();
    % d = 0;
    for i = 1:5
        subd(i) = detv(B(2:end,cp(i,:)));
        d       = (-1)^(rr(i)-1) * B(1,rr(i))*subd(i) + d;
    end
end    