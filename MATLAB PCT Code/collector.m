clc; clear; clf;
for j=1:9;
     A(:,:,j)=importdata([num2str(j) '.dat']);
end
n=1; L=length(A);
for j=[1 4 7]
    B(:,:,n)=horzcat(A(2:L-1,2:L-1,j),A(2:L-1,2:L-1,j+1),A(2:L-1,2:L-1,j+2));
    n=n+1;
end
C=vertcat(B(:,:,1),B(:,:,2),B(:,:,3));
image_hex(C)
axis equal
axis tight
