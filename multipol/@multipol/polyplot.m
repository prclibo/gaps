function polyplot(polynomial)
%plots the monomials of the polynomial
%only available for polynomials of up to three variables
monoms = polynomial.monomials;
if(size(monoms, 1) > 3); error('to many variables');end;

if(size(monoms, 1) == 3)
    N = max(monoms(end,:));
    for k = 0 : N
        submonoms = monoms(:, find(monoms(end,:)==k));
        subplot(1,N+1,k+1);
        plot(submonoms(1,:),submonoms(2,:),'*r');
    end
else
    if(size(monoms, 1) == 2)
        plot(monoms(1,:),submonoms(2,:),'*r');
    else
        plot(monoms, '*r');
    end
end
