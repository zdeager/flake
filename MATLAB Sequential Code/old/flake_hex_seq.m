% flake_seq_hex.m
% Sequential MATLAB implementation of Reiter's local cellular model for snow crystal growth
%
% zde

clear;clc;clf
L=10;T=3; % matrix size, # of steps
alpha=1;beta=0.35;gamma=0.001; %parameters

%STEP1:initialize matrices and put seed in center
A=beta*ones(L);A2=A;A2n=A2;A1=zeros(L);
c=ceil(L/2);A(c,c)=alpha;

for n=1:T
    for i=2:L-1 %STEP2:scan cells
        for j=2:L-1
            if mod(i,2) == 1 % if row is odd, neighbors are (i-1,j-1),(i-1,j),(i,j+1),(i+1,j),(i+1,j-1),(i,j-1)
                if A(i,j)>=alpha||A(i-1,j-1)>=alpha||A(i-1,j)>=alpha||A(i,j+1)>=alpha...
            	||A(i+1,j)>=alpha||A(i+1,j-1)>=alpha||A(i,j-1)>=alpha
                    A1(i,j)=A(i,j)+gamma; %STEP3:grow ice
                    A2(i,j)=0;
            	else
                    A1(i,j)=0; 
                    A2(i,j)=A(i,j);
                end
            else % else row is even, neighbors are (i-1,j+1),(i,j+1),(i+1,j+1),(i+1,j),(i,j-1),(i-1,j)
            	if A(i,j)>=alpha||A(i-1,j+1)>=alpha||A(i,j+1)>=alpha||A(i+1,j+1)>=alpha...
            	||A(i+1,j)>=alpha||A(i,j-1)>=alpha||A(i-1,j)>=alpha
                    A1(i,j)=A(i,j)+gamma; %STEP3:grow ice
                    A2(i,j)=0;
            	else
                    A1(i,j)=0; 
                    A2(i,j)=A(i,j);
                end
            end
        end
    end

    for i=2:L-1 %STEP4:diffuse water
        for j=2:L-1
            if mod(i,2) == 1 % if row is odd, neighbors are (i-1,j-1),(i-1,j),(i,j+1),(i+1,j),(i+1,j-1),(i,j-1)
                avg_neigh = (A2(i-1,j-1) + A2(i-1,j) + A2(i,j+1) + A2(i+1,j) + A2(i+1,j-1) + A2(i,j-1)) / 6;
            else % else row is even, neighbors are (i-1,j+1),(i,j+1),(i+1,j+1),(i+1,j),(i,j-1),(i-1,j)
                avg_neigh = (A2(i-1,j+1) + A2(i,j+1) + A2(i+1,j+1) + A2(i+1,j) + A2(i,j-1) + A2(i-1,j)) / 6;
            end
            A2n(i,j) = (A2(i,j) + avg_neigh) / 2; % weighted average of current cell and 6 neighbors
        end
    end

    %image_hex(A);drawnow; pause(1); %plot hex growth
    A=A1+A2n; %STEP5:add updated water and ice
    A2=A2n; %update water for next iteration
end
