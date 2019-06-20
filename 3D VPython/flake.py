from vpython import *

W = 60 # width
H = 60 # height
D = 60 # depth
# determine number of neighbors
if D > 1:
  N = 12
else:
  N = 6
R = 1 # radius

T = 15 # num of iterations

# growth params
alpha = 1
beta  = 0.9
gamma = 0.05

# cell color (lt. blue-ish)
color = vector(204/255, 1, 1)

def get_neighbors(idx):
  # see https://stackoverflow.com/questions/9982458/creating-a-sphere-packing-nearest-neighbor-list-from-integers
  A = W * H

  plane = int(idx / A)
  plane_index = idx % A
  row = int(plane_index / W)
  col = plane_index % W

  r = -1 if row % 2 else 1   # (-1)**row
  p = -1 if plane % 2 else 1 # (-1)**plane

  nbors = []

  # first include neighbors in same plane
  if col != W-1: nbors.append(idx+1)
  if col != 0:   nbors.append(idx-1)
  if row != H-1: nbors.append(idx+W)
  if row != 0:   nbors.append(idx-W)
  if (col != 0 or r > 0) and (col != W-1 or r < 0):
    if row != H-1: nbors.append(idx+W+r)
    if row != 0:   nbors.append(idx-W+r)

  # now add neighbors from other planes
  if plane != D-1: nbors.append(idx+A)
  if plane != 0:   nbors.append(idx-A)

  if (col != 0 or p < 0) and (col != W-1 or p > 0):
    if plane != D-1: nbors.append(idx+A-p)
    if plane != 0:   nbors.append(idx-A-p)

  if ((col != W - 1 or p > 0 or r < 0) and
      (col != 0 or p < 0 or r > 0) and
      (row != H-1 or p < 0) and
      (row != 0 or p > 0)):
    if plane != D-1:
      nbors.append(idx + A + p*W + int((r-p)/2)) #10
    if plane != 0:
      nbors.append(idx - A + p*W + int((r-p)/2)) #11

  return nbors

def get_pos(idx):
  k = int(idx / (W*H));
  idx -= (k*W*H);
  j = int(idx / W);
  i = idx % W;
  #print(i,j,k)
  x = 2*i + ((j+k) % 2)
  y = sqrt(3)*(j + (1/3)*(k % 2))
  z = ((2*sqrt(6))/3)*k

  return vector(x,y,z)*R

cells = []

# initialize 
for idx in range(W*H*D):
  cell = sphere()
  cell.shininess = 0
  cell.pos = get_pos(idx)
  cell.radius = R
  cell.A = beta
  cell.A2 = cell.A
  cell.A2n = cell.A2
  cell.A1 = 0
  cell.neighbors = get_neighbors(idx)
  cells.append(cell)
  cell.visible = False
center_idx = int(D/2)*(H*W) + int(H/2)*W + int(W/2)
cells[center_idx].A = alpha

# iterate
for iter in range(T):
  print(iter)
  # grow ice
  for idx in range(W*H*D):
    neighbors = cells[idx].neighbors
    if len(neighbors) != N: 
      continue
    # determine receptive
    receptive = False 
    for neighbor in neighbors:
      if cells[neighbor].A >= alpha or cells[idx].A >= alpha:
        receptive = True
        break
    if receptive:
      cells[idx].A1 = cells[idx].A + gamma
      cells[idx].A2 = 0
    else:
      cells[idx].A1 = 0
      cells[idx].A2 = cells[idx].A 

  # diffuse
  for idx in range(W*H*D):
    neighbors = cells[idx].neighbors
    if len(neighbors) != N: 
      continue
    avg = 0
    for neighbor in neighbors:
      avg += cells[neighbor].A2
    avg /= N 
    cells[idx].A2n = (cells[idx].A2 + avg) / 2

  for idx in range(W*H*D):
    cells[idx].A = cells[idx].A1 + cells[idx].A2n #add updated water and ice
    cells[idx].A2 = cells[idx].A2n # update water for next step

# color cells 
#scene.userzoom = False
#scene.userspin = False
scene.center = get_pos(center_idx)
scene.width = 600
scene.height = 600

A_max = -float("inf")
A_min = float("inf")
for idx in range(W*H*D):
  if cells[idx].A > A_max:
    A_max = cells[idx].A
  if cells[idx].A < A_min:
    A_min = cells[idx].A
for idx in range(W*H*D):
  norm = (cells[idx].A-A_min) * (1 / (A_max-A_min));
  cells[idx].color = color*norm
  cells[idx].opacity = norm-.1
  if norm > .5:
    cells[idx].visible = True

