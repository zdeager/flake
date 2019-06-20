% image_hex.m
% Draw the hexagonal grid represented by the matrix A
% NOTE: this implementation computes the vertices of each hexagon and makes a single call to patch - FAST
% zde

function [] = image_hex(A)
	clf;clc
	colormap(bone(20));
	colorbar
	th = linspace(0,2*pi,6+1).'+pi/6; th = th(2:end);
	xh = 1/sqrt(3)*cos(th);
	yh = 1/sqrt(3)*sin(th);

	N = size(A,1);
	xx = repmat(1:N,6,1) + repmat(xh,[1 N]);
	yy = repmat(yh,[1 N]);

	tmp=repmat(reshape(repmat(1:2,N,1),1,[]),1,floor(N/2));
    if mod(size(A),2) == 1
        tmp = [tmp tmp(:,1:N)];
    end
	xx = .5*repmat(tmp,6,1) + repmat(xx,[1 N]);
	yy = sqrt(3)/2*repmat(reshape(repmat(1:N,N,1),1,[]),6,1) + repmat(yy,[1 N]);
    
	tmp = A';
    patch(xx,yy,tmp(:), 'edgecolor','none');
    daspect([1 1 1]);