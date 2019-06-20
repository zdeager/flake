% flake_seq.m
% Sequential MATLAB implementation of modified Reiter's local 
% cellular model for snow crystal growth with square cells instead of hexagonal cells
%
% zde

clear;clc;clf;colormap bone
L=15;T=3;CI=60; %grid size,# of steps, and color map intesnsity
alpha=.5;beta=0.2;gamma=0.0025; %parameters
%STEP1:initialize grid
A=beta*ones(L);A2=A;A2n=A2;A1=zeros(L);
c=ceil(L/2);A(c,c)=alpha;
for n=1:T
    for i=2:L-1%STEP2:scan cells
        for j=2:L-1
            if A(i,j)>=alpha||A(i+1,j)>=alpha||A(i-1,j)>=alpha||A(i,j+1)>=alpha...
            ||A(i,j-1)>=alpha||A(i+1,j+1)>=alpha||A(i-1,j-1)>=alpha||A(i-1,j+1)>=alpha...
            ||A(i+1,j-1)>=alpha
                A1(i,j)=A(i,j)+gamma;%STEP3:grow ice
                A2(i,j)=0;
            else
                A1(i,j)=0;
                A2(i,j)=A(i,j);
            end
        end
    end
    for i=2:L-1%STEP4:diffuse water
        for j=2:L-1
            A2n(i,j)=(A2(i+1,j)+A2(i-1,j)+A2(i,j+1)+A2(i,j-1)...
            +A2(i+1,j+1)+A2(i-1,j-1)+A2(i+1,j-1)+A2(i-1,j+1))/8;
        end
    end
    %image(CI*A);drawnow; pause(.5); %plot growth
    A=A1+A2n;%STEP5:add updated water and ice
    A2=A2n;%update water for next iteration
end
