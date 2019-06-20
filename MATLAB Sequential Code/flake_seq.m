% faster version of flake_seq using convolutions
% zde
clc;clear;

L=12;
T=1; % matrix size, # of steps
alpha=1.00001;
beta=0.95;
gamma=0.001; %parameters

A=beta*ones(L);
c=ceil(L/2);A(c,c)=alpha;

tic
for tt=1:T
    tmp = conv2(1.0*(A>=alpha),[1 1 1; 1 1 1; 1 1 1],'same');
    A1 = (1.0*(tmp>=1)).*(A+gamma);
    A2 = (1.0*(tmp==0)).*A;
    
    A2avg = conv2(A2,[1 1 1; 1 0 1; 1 1 1],'same')/8;
    A2n = (A2+A2avg)/2;
    
    A = A1 + A2n;
    A2 = A2n;
end
toc
    
    
    
