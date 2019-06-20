% flake_hex_MPI_sq.m
% Parallel MATLAB implementation of modified Reiter's local cellular model for snow crystal growth
% USING MATLAB'S PCT
%
% nsb, zde

% Square based domain decomp.
% L : dimension of matrix
% T : num of iterations
% alpha,beta,gamma : inital parameters 
function time=flake_hex_MPI_sq(L, T, alpha, beta, gamma) 
    matlabpool 9;
    spmd

    rt = sqrt(numlabs); % number of rows(cols) of nodes
    cols=ceil(L/rt)+2; rows=cols; % num of rows/cols for matrix on this node
    A=beta*ones(rows,cols); A1=0*A;A2=A;A2n=A2; % initalize matrix on this node
    
    c = ceil(L/2); % compute entry for seed
    col = ceil(c/(rows-2)); row = ceil(c/(cols-2)); % compute col/row of node seed belongs in
    seed = ((row - 1)*rt)+col; % node seed belongs in
    if labindex == seed
        if mod(numlabs,2) == 1 % number of nodes is odd
            A(ceil(rows/2), ceil(cols/2)) = alpha; % place seed in center of this node
        else
            A(rows-1,cols-1) = alpha; % place seed in bottom-left of this node
        end
    end; tic	
    for n=1:T
        %% send updated A matrix between partitions
        if labindex == 1 % top-left node
            labSend(A(:,cols-1),labindex+1); % send last col to right node
            labSend(A(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            labSend(A(rows-1,:),labindex+rt); % send bottom row to bottom node
            A(:,cols)=labReceive(labindex+1); % get first col from right node
            A(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
        elseif labindex == rt % top-right node
            labSend(A(:,2),labindex-1); % send first col to left node
            labSend(A(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            labSend(A(rows-1,:),labindex+rt); % send bottom row to bottom node
            A(:,1)=labReceive(labindex-1); % get last col from left node
            A(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
        elseif labindex == numlabs - (rt - 1) % bottom-left node
            labSend(A(:,cols-1),labindex+1); % send last col to right node
            labSend(A(2,cols-1),labindex-rt+1); % send top-right element to top-right node 
            labSend(A(2,:),labindex-rt); % send top row to top node
            A(:,cols)=labReceive(labindex+1); % get first col from right node
            A(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
        elseif labindex == numlabs % bottom-right node
            labSend(A(:,2),labindex-1); % send first col to left node
            labSend(A(2,2),labindex-rt-1); % send top-left element to top-left node 
            labSend(A(2,:),labindex-rt); % send top row to top node
            A(:,1)=labReceive(labindex-1); % get last col from left node
            A(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
        elseif labindex > 1 && labindex < rt % top node
            labSend(A(:,cols-1),labindex+1); % send last col to right node
            labSend(A(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            labSend(A(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A(:,2),labindex-1); % send first col to left node
            labSend(A(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            A(:,cols)=labReceive(labindex+1); % get first col from right node
            A(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A(:,1)=labReceive(labindex-1); % get last col from left node
            A(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
            A(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
        elseif mod(labindex, rt) == 1 % left node
            labSend(A(2,:),labindex-rt); % send top row to top node
            labSend(A(:,cols-1),labindex+1); % send last col to right node
            labSend(A(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A(2,cols-1),labindex-rt+1); % send top-right element to top-right node
            labSend(A(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            A(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A(:,cols)=labReceive(labindex+1); % get first col from right node
            A(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
            A(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
        elseif mod(labindex, rt) == 0 % right node
            labSend(A(:,2),labindex-1); % send first col to left node
            labSend(A(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A(2,:),labindex-rt); % send top row to top node
            labSend(A(2,2),labindex-rt-1); % send top-left element to top-left node
            labSend(A(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            A(:,1)=labReceive(labindex-1); % get last col from left node
            A(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
            A(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
        elseif labindex > numlabs - (rt - 1) && labindex < numlabs % bottom node
        	labSend(A(2,:),labindex-rt); % send top row to top node
            labSend(A(:,2),labindex-1); % send first col to left node
            labSend(A(:,cols-1),labindex+1); % send last col to right node
            labSend(A(2,cols-1),labindex-rt+1); % send top-right element to top-right node 
            labSend(A(2,2),labindex-rt-1); % send top-left element to top-left node
            A(:,1)=labReceive(labindex-1); % get last col from left node
            A(:,cols)=labReceive(labindex+1); % get first col from right node
            A(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
            A(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
        else % middle node
        	labSend(A(2,:),labindex-rt); % send top row to top node
            labSend(A(:,2),labindex-1); % send first col to left node
            labSend(A(:,cols-1),labindex+1); % send last col to right node
            labSend(A(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A(2,cols-1),labindex-rt+1); % send top-right element to top-right node 
            labSend(A(2,2),labindex-rt-1); % send top-left element to top-left node
            labSend(A(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            labSend(A(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            A(:,1)=labReceive(labindex-1); % get last col from left node
            A(:,cols)=labReceive(labindex+1); % get first col from right node
            A(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
            A(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
            A(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
            A(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
        end

	offset = (floor((labindex - 1)/rt) * (rows - 2)) - 1; % compute offset for row in agregate matrix

        for i=2:rows-1 %STEP 2: scan cells
            for j=2:cols-1
                if mod(i+offset,2) == 1 % if row is odd, neighbors are (i-1,j-1),(i-1,j),(i,j+1),(i+1,j),(i+1,j-1),(i,j-1)
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
        
        %% send updated A2 matrix between partitions
        if labindex == 1 % top-left node
            labSend(A2(:,cols-1),labindex+1); % send last col to right node
            labSend(A2(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            labSend(A2(rows-1,:),labindex+rt); % send bottom row to bottom node
            A2(:,cols)=labReceive(labindex+1); % get first col from right node
            A2(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A2(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
        elseif labindex == rt % top-right node
            labSend(A2(:,2),labindex-1); % send first col to left node
            labSend(A2(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            labSend(A2(rows-1,:),labindex+rt); % send bottom row to bottom node
            A2(:,1)=labReceive(labindex-1); % get last col from left node
            A2(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A2(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
        elseif labindex == numlabs - (rt - 1) % bottom-left node
            labSend(A2(:,cols-1),labindex+1); % send last col to right node
            labSend(A2(2,cols-1),labindex-rt+1); % send top-right element to top-right node 
            labSend(A2(2,:),labindex-rt); % send top row to top node
            A2(:,cols)=labReceive(labindex+1); % get first col from right node
            A2(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A2(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
        elseif labindex == numlabs % bottom-right node
            labSend(A2(:,2),labindex-1); % send first col to left node
            labSend(A2(2,2),labindex-rt-1); % send top-left element to top-left node 
            labSend(A2(2,:),labindex-rt); % send top row to top node
            A2(:,1)=labReceive(labindex-1); % get last col from left node
            A2(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A2(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
        elseif labindex > 1 && labindex < rt % top node
            labSend(A2(:,cols-1),labindex+1); % send last col to right node
            labSend(A2(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            labSend(A2(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A2(:,2),labindex-1); % send first col to left node
            labSend(A2(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            A2(:,cols)=labReceive(labindex+1); % get first col from right node
            A2(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A2(:,1)=labReceive(labindex-1); % get last col from left node
            A2(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
            A2(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
        elseif mod(labindex, rt) == 1 % left node
            labSend(A2(2,:),labindex-rt); % send top row to top node
            labSend(A2(:,cols-1),labindex+1); % send last col to right node
            labSend(A2(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A2(2,cols-1),labindex-rt+1); % send top-right element to top-right node
            labSend(A2(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            A2(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A2(:,cols)=labReceive(labindex+1); % get first col from right node
            A2(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A2(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
            A2(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
        elseif mod(labindex, rt) == 0 % right node
            labSend(A2(:,2),labindex-1); % send first col to left node
            labSend(A2(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A2(2,:),labindex-rt); % send top row to top node
            labSend(A2(2,2),labindex-rt-1); % send top-left element to top-left node
            labSend(A2(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            A2(:,1)=labReceive(labindex-1); % get last col from left node
            A2(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A2(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A2(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
            A2(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
        elseif labindex > numlabs - (rt - 1) && labindex < numlabs % bottom node
        	labSend(A2(2,:),labindex-rt); % send top row to top node
            labSend(A2(:,2),labindex-1); % send first col to left node
            labSend(A2(:,cols-1),labindex+1); % send last col to right node
            labSend(A2(2,cols-1),labindex-rt+1); % send top-right element to top-right node 
            labSend(A2(2,2),labindex-rt-1); % send top-left element to top-left node
            A2(:,1)=labReceive(labindex-1); % get last col from left node
            A2(:,cols)=labReceive(labindex+1); % get first col from right node
            A2(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A2(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
            A2(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
        else % middle node
        	labSend(A2(2,:),labindex-rt); % send top row to top node
            labSend(A2(:,2),labindex-1); % send first col to left node
            labSend(A2(:,cols-1),labindex+1); % send last col to right node
            labSend(A2(rows-1,:),labindex+rt); % send bottom row to bottom node
            labSend(A2(2,cols-1),labindex-rt+1); % send top-right element to top-right node 
            labSend(A2(2,2),labindex-rt-1); % send top-left element to top-left node
            labSend(A2(rows-1,cols-1),labindex+rt+1); % send bottom-right element to bottom-right node 
            labSend(A2(rows-1,2),labindex+rt-1); % send bottom-left element to bottom-left node 
            A2(:,1)=labReceive(labindex-1); % get last col from left node
            A2(:,cols)=labReceive(labindex+1); % get first col from right node
            A2(1,:)=labReceive(labindex-rt); % get bottom row from top node
            A2(rows,:)=labReceive(labindex+rt); % get top row from bottom node
            A2(1,cols)=labReceive(labindex-rt+1); % get bottom-left element from top-right node
            A2(1,1)=labReceive(labindex-rt-1); % get bottom-right element from top-left node
            A2(rows,cols)=labReceive(labindex+rt+1); % get top-left element from bottom-right node
            A2(rows,1)=labReceive(labindex+rt-1); % get top-right element from bottom-left node
        end
        
        for i=2:rows-1 %STEP 4: diffuse water
            for j=2:cols-1
                if mod(i + offset,2) == 1 % if row is odd, neighbors are (i-1,j-1),(i-1,j),(i,j+1),(i+1,j),(i+1,j-1),(i,j-1)
                    avg_neigh = (A2(i-1,j-1) + A2(i-1,j) + A2(i,j+1) + A2(i+1,j) + A2(i+1,j-1) + A2(i,j-1)) / 6;
                else % else row is even, neighbors are (i-1,j+1),(i,j+1),(i+1,j+1),(i+1,j),(i,j-1),(i-1,j)
                    avg_neigh = (A2(i-1,j+1) + A2(i,j+1) + A2(i+1,j+1) + A2(i+1,j) + A2(i,j-1) + A2(i-1,j)) / 6;
                end
                A2n(i,j) = (A2(i,j) + avg_neigh) / 2; % weighted average of current cell and 6 neighbors
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
end
