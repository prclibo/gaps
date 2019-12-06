function B=adj(A)

n=size(A,1); d=1:n-1;
B=A(1,1)*zeros(n); AA=[A,A;A,A]';

if mod(n,2)==1,
    for j=1:n
        for k=1:n
            B(j,k)=det(AA(j+d,k+d));
        end
    end
else
    for j=1:n
        for k=1:n
            B(j,k)=(-1)^(j+k)*det(AA(j+d,k+d));
        end
    end
end
