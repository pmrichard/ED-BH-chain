function [Cf] = condFrac(G)
%Compute the condensate fraction.

global D;
global basis;
global index;
global M;
global N;
global factor

dens = zeros(M, M);
for m = 1 : D
    if G(m) ~= 0
        dens(1,1) = dens(1,1)+G(m)*G(m)*basis(1,m);
        for s1 = 2 : M
            if (basis(1,m) > 0)
                final = basis(:,m);
                final(s1) = final(s1)+1;
                final(1) = final(1)-1;
                value = factor*final;
                finalnum = index(search(value));
                if G(finalnum) ~= 0;
                    dens(s1,1) = dens(s1,1)+G(m)*G(finalnum)*sqrt(final(s1)*(final(1)+1));
                end
            end            
        end
    end
end

for m = 1 : M-1
    for k = 1 : M-m
        dens(m+k,1+k) = dens(m,1);
    end
end
for m = 2 : M
    for k = 1 : m-1
        dens(k,m) = dens(m,k);
    end
end
    
Cf = eigs(dens,1,'la')/N;

end