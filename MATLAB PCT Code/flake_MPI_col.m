% flake_MPI.m
% Parallel MATLAB implementation of modified Reiter's local 
% cellular model for snow crystal growth with square cells instead of hexagonal cells
% USING MATLAB'S PCT
%
% nsb, zde

% Column based domain decomp.
% L : dimension of matrix
% T : num of iterations
% alpha,beta,gamma : inital parameters 
function time=flakeMPI(L, T, alpha, beta, gamma) 
    cols=ceil(L/numlabs)+2; rows=ceil(L/numlabs)*numlabs;
    A=beta*ones(rows,cols); A1=0*A;A2=A;A2n=A2; 
    middle1=floor(median(1:numlabs)); middle2=ceil(median(1:numlabs));
    if middle1~=middle2 %STEP 1: initialize grid
        if labindex==middle1, A(ceil(rows/2),cols-1)=alpha; end
    else
        if labindex==middle1, A(ceil(rows/2),ceil(cols/2))=alpha; end
    end; tic
    for n=1:T
        %% send updated A matrix between partitions
        if labindex>1, labSend(A(:,2),labindex-1); end
        if labindex<numlabs
            labSend(A(:,cols-1),labindex+1);
            A(:,cols)=labReceive(labindex+1); 
        end
        if labindex>1, A(:,1)=labReceive(labindex-1); end
	    for i=2:rows-1 %STEP 2: scan cells
            for j=2:cols-1
                if A(i,j)>=alpha||A(i+1,j)>=alpha||A(i-1,j)>=alpha||A(i,j+1)>=alpha...
                        ||A(i,j-1)>=alpha||A(i+1,j+1)>=alpha||A(i-1,j-1)>=alpha||A(i-1,j+1)>=alpha...
                        ||A(i+1,j-1)>=alpha
                    A1(i,j)=A(i,j)+gamma; A2(i,j)=0; %STEP 3: grow ice
                else
                    A1(i,j)=0; A2(i,j)=A(i,j);
                end
            end
        end
        %% send updated A2 matrix between partitions
        if labindex>1, labSend(A2(:,2),labindex-1); end
        if labindex<numlabs
            labSend(A2(:,cols-1),labindex+1);
            A2(:,cols)=labReceive(labindex+1); 
        end
        if labindex>1, A2(:,1)=labReceive(labindex-1); end
        for i=2:rows-1 %STEP 4: diffuse water
            for j=2:cols-1
                A2n(i,j)=(A2(i+1,j)+A2(i-1,j)+A2(i,j+1)+A2(i,j-1)...
                    +A2(i+1,j+1)+A2(i-1,j-1)+A2(i+1,j-1)+A2(i-1,j+1))/8;
            end
        end
        A=A1+A2n; A2=A2n; %combine and update for next iteration
	
    %   filename=strcat('it',int2str(n),'.jpg'); 
    %   AA=gop(@horzcat,A(:,2:cols-1));  
    %   if labindex==1
    %       imwrite(AA/2,filename);
    %   end
    end;
    %if labindex==1
    %    imwrite(AA/2,'flake.jpg')
    %end
    %imwrite(A,[num2str(labindex),'.jpg']);
    dlmwrite([num2str(labindex),'.dat'], A);
    
    time=toc;
end
