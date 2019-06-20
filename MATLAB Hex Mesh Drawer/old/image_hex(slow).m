% image_hex.m
% Draw the hexagonal grid represented by the matrix A
% NOTE: this is the first implementation that plots each patch iteratively - SLOW
% zde
function [] = image_hex(A)

	clf;clc
	
	% vertices for hexagon inscribed in a circle with radius 
	% r = 50, centered at origin (rotated 90 deg.)
	r = 1;
	x0 = [0,r*-sqrt(3)/2,r*-sqrt(3)/2,0,r*sqrt(3)/2,r*sqrt(3)/2]; x = x0;
	y0 = [r*-1,r*-.5,r*.5,r*1,r*.5,r*-.5]; y = y0;

	normA = mat2gray(A); % normalize the entries of A for color intensity

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
			patch(x,y,'EdgeColor','b','FaceColor',normA(i,j)*[1,1,1]) % draw hex grid		
		end
	end

