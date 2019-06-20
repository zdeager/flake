% image_hex.m
% Draw the hexagonal grid represented by the matrix A
% NOTE: this implementation computes the vertices and faces of each hexagon and makes a single call to patch - FASTER
% zde

function [] = image_hex(A)
	tic
	clf;clc
	
	verts = zeros(length(A)^2 * 6,2); % store x,y-coords of each hexagons' vertices
	faces = zeros(length(A)^2, 6); % stores which vertices corespond to each hexagon
	color_mat = zeros(length(A)^2,3); % stores rgb color of each hexagon

	% vertices for hexagon inscribed in a circle with radius r centered at origin (rotated 90 deg.)
	r = 1;
	x0 = [0,r*-sqrt(3)/2,r*-sqrt(3)/2,0,r*sqrt(3)/2,r*sqrt(3)/2]; x = x0;
	y0 = [r*-1,r*-.5,r*.5,r*1,r*.5,r*-.5]; y = y0;

	normA = mat2gray(A); % normalize the entries of A for color intensity

	% compute coords for each hexagons' vertices
	k = 1; % looping variable for hexagons' vertices and color matrices
	for i = 1:length(A)
		for j = 1:length(A)
			if j == 1 && mod(i,2) == 1 % starting a new odd row, shift down
				x = x0; % reset x vertices
				y -= r * 3/2; % shift down by 3/2 * r
			elseif j == 1 && mod(i,2) == 0 % if starting a new even row, shift over and down
				x = x0; % reset x vertices
				x += r * sqrt(3)/2; % shift over by sqrt(3)/2 * r
				y -= r * 3/2; % shift down by 3/2 * r
			else % shift over
				x += r * sqrt(3); % shift over by sqrt(3)*r
			end
			verts((k-1)*6+1:(k-1)*6+6,1) = x;
			verts((k-1)*6+1:(k-1)*6+6,2) = y;
			faces(k,:) = [(k-1)*6+1:(k-1)*6+6];
			color_mat(k,:) = normA(i,j)*[1,1,1]; % compute hexagon's fill color
			k = k + 1;
		end
	end
	patch('Faces', faces, 'Vertices', verts, 'EdgeColor','none','FaceVertexCData',color_mat); % draw patches
	toc

