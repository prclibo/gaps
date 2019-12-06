function roots = zp_sqrt(n, p)

% https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python

if legendre(n, p) ~= 1
    roots = [];
    return;
else
    q = p - 1;
    s = 0;
    while logical(mod(q, 2) == 0)
        q = q / 2;
        s = s + 1;
    end
    if logical(s == 1)
        roots = powermod(n, (p + 1) / 4, p);
        return;
    else
        for z = 2:p - 1
            if p - 1 == legendre(z, p), break; end
        end
        c = powermod(z, q, p);
        r = powermod(n, (q + 1) / 2, p);
        t = powermod(n, q, p);
        m = s;
        t2 = 0;
        while logical(mod((t - 1), p) ~= 0)
            t2 = powermod(t, 2, p);
            for i = 1:m
                if logical(mod((t2 - 1), p) == 0), break; end
                t2 = powermod(t2, 2, p);
            end
            b = powermod(c, 2^(m - i - 1), p);
            r = mod(r * b, p);
            c = mod(b * b, p);
            t = mod(t * c, p);
            m = i;
        end
        roots = [r, p - r];
        return
    end
end
end
function ret = legendre(a, p)
ret = powermod(a, (p - 1) / 2, p);
end