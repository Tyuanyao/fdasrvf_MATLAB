function basis = findBasisNormal(q)
% Return basis vectors for normal space at q. basis is a cell array of size

[n,T]=size(q);
fs = zeros(n,T,n);

for i = 1:T
    for j = 1:n
        e = zeros(n,1);
        e(j) = 1;
        fs(:,i,j) = q(1,i)*q(:,i)/norm(q(:,i))+norm(q(:,i))*e;

    end
end

integrandb = zeros(T,n);

for i=1:T
    for j = 1:n
        integrandb(i,j) = q(:,i)'*fs(:,i,j);
        integrandb(i,j) = q(:,i)'*fs(:,i,j);
    end
end

basis = cell(3, 1);
for j = 1:n
    basis{j} = fs(:,:,j) - q*trapz(linspace(0,1,T),integrandb(:,j));
end

