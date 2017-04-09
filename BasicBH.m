%The simplest BH model of M sites, N is the number of bosons.
clear; close all; 

global D; %Dimension of the Hilbert space.
global basis; %Basis vectors
global address; %Hashing address
global index; %Oringin position
global M;
global N;
global factor; %Hash factor

Min = 5;
Max = 8;
for M = Min : Max
    tic
N = M;
D = factorial(N+M-1)/(factorial(N)*factorial(M-1));
address = zeros(1,D);

basis = zeros(M,D); 
basis(1,1) = N;

factor = sqrt(100*(1:M)+3); 
address(1) = factor*basis(:,1);
s = 1; %Update the Hashing address as well as basis.
while (basis(M,s) < N)
    s1 = M-1;
    while (basis(s1,s)==0)
        s1 = s1-1;
    end
    for s2 = 1 : s1-1
        basis(s2,s+1) = basis(s2,s);
    end
    basis(s1,s+1) = basis(s1,s)-1;
    basis(s1+1,s+1) = N-sum(basis(1:s1,s+1));
    address(s+1) = factor*basis(:,s+1);     
    s = s+1;
end
[address,index] = sort(address); %Quicksort the addresses.

%Generate the Harmiltonian
gendiag = zeros(2,D);
for s1 = 1 : D
    gendiag(1,s1) = s1;
    gendiag(2,s1) = 0.5*sum(basis(:,s1).^2-basis(:,s1));
end
Int = sparse(gendiag(1,:),gendiag(1,:),gendiag(2,:),D,D); %H_int

genoff = zeros(3,D*M*2);
s = 0;
target=[(2:M),1];
for s1 = 1 : D   
    for s2 = 1 : M
         if (basis(s2,s1) > 0)
             s = s+1;
             final = basis(:,s1);
             final(s2) = final(s2)-1;
             s3 = target(s2);
             final(s3) = final(s3)+1;
             value = factor*final;
             addindex = search(value);
             
             genoff(1,s) = s1;
             genoff(2,s) = index(addindex);
             genoff(3,s) = -sqrt(final(s3)*(final(s2)+1));
         end
    end
end
genoff = genoff(:,1:s);
Kin = sparse(genoff(1,:),genoff(2,:),genoff(3,:),D,D);
Kin = Kin + Kin'; %H_kin

%Measure
Vlist=0:0.5:20;
Cf = zeros(1,length(Vlist));
for s1 = 1:length(Vlist)
    H = Kin + Vlist(s1)*Int;
    [V,d] = eigs(H, 1, 'sa'); %The ground state
    Cf(s1) = condFrac(V);    
end

%Plot
plot(Vlist,Cf);
ylim([0 1])
xlim([0 max(Vlist)])
xlabel('V/J','fontsize',14)
ylabel('f_c','fontsize',14)
title('Condensate fraction')
hold on

M
toc
end
legend(num2str((Min:Max)','N = %d'));
