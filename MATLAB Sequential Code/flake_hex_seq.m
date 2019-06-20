% faster version of flake_hex_seq using convolutions
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
    tmp_even = conv2(1.0*(A>=alpha),[0 1 1; 1 1 1; 0 1 1],'same');
    tmp_odd = conv2(1.0*(A>=alpha),[1 1 0; 1 1 1; 1 1 0],'same');
    A1 = zeros(size(A)); A2 = zeros(size(A));
    tmp = tmp_odd>=1; A1(2:2:end,:) = tmp(2:2:end,:).*(A(2:2:end,:)+gamma);
    tmp = tmp_odd==0; A2(2:2:end,:) = tmp(2:2:end,:).*A(2:2:end,:);
    tmp = tmp_even>=1; A1(1:2:end,:) = tmp(1:2:end,:).*(A(1:2:end,:)+gamma);
    tmp = tmp_even==0; A2(1:2:end,:) = tmp(1:2:end,:).*A(1:2:end,:);
    
    tmp_even = conv2(A2,[0 1 1; 1 0 1; 0 1 1],'same')/6;
    tmp_odd = conv2(A2,[1 1 0; 1 0 1; 1 1 0],'same')/6;
    A2avg = zeros(size(A2));
    A2avg(2:2:end,:) = tmp_odd(2:2:end,:);
    A2avg(1:2:end,:) = tmp_even(1:2:end,:);
    A2n = (A2+A2avg)/2;
    
    A = A1 + A2n;
    A2 = A2n;
end
toc