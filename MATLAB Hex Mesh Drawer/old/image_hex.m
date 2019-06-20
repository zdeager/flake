% image_hex.m
% Draw the hexagonal grid represented by the matrix A
% NOTE: this implementation computes the vertices of each hexagon and makes a single call to patch - FAST
% zde

function [] = image_hex(A)
	clf;clc
	rows = length(A(:,1))
	cols = length(A(1,:))
	patches_x = zeros(6,rows*cols); % stores x-coords of hexagons' vertices
	patches_y = zeros(6,rows*cols); % stores y-coords of hexagons' vertices
	color_mat = zeros(1,rows*cols,3); % stores rgb color of each hexagon's face

	% vertices for hexagon inscribed in a circle with radius r centered at origin (rotated 90 deg.)
	r = 1;
	x0 = [0,r*-sqrt(3)/2,r*-sqrt(3)/2,0,r*sqrt(3)/2,r*sqrt(3)/2]; x = x0;
	y0 = [r*-1,r*-.5,r*.5,r*1,r*.5,r*-.5]; y = y0;

	normA = mat2gray(A); % normalize the entries of A for color intensity

	% compute coords for each hexagons' vertices
	k = 1; % looping variable for hexagons' vertices and color matrices
	for i = 1:rows
		for j = 1:cols
			if j == 1 && mod(i,2) == 1 % starting a new odd row, shift down
				x = x0; % reset x vertices
				y = y - r * 3/2; % shift down by 3/2 * r
			elseif j == 1 && mod(i,2) == 0 % if starting a new even row, shift over and down
				x = x0; % reset x vertices
				x = x + r * sqrt(3)/2; % shift over by sqrt(3)/2 * r
				y = y - r * 3/2; % shift down by 3/2 * r
			else % shift over
				x = x + r * sqrt(3); % shift over by sqrt(3)*r
			end
			color_mat(1,k,:) = normA(i,j)*[.8 1 1]; % compute hexagon's fill color
			patches_x(:,k) = x; % add x coords of hexagon's vertices to matrix
			patches_y(:,k) = y; % add y coords of hexagon's vertices to matrix
			k = k + 1;
		end
	end
	patch(patches_x, patches_y, color_mat, 'EdgeColor', 'none'); % draw patches
