/**
 * A Processing implementation Reiter's cellular model for snow crystal growth
 * flake.pde
 * zde
 *
 * Press SPACE BAR to pause/play simulation
 * Press R to reset the cells and pause
 * On reset picks new values for beta/gamma at random
 * On reset/pause, click a cell to place a seed
 *
 * NOTE: currently using square lattice
 * TODO: - use hexagonal lattice
 */

// Size of cells (px)
final int cellSize = 5;

// Simulation parameters
float alpha;
float beta;
float gamma;

// Fill color (lt. blue-ish)
final float R = 204;
final float G = 255;
final float B = 255;

// Variables for timer
int interval = 1000;
int lastRecordedTime = 0;

// Matrices
float[][] A; 
float[][] A1; 
float[][] A2; 
float[][] A2n;
float A_max,A_min;

// Pause
boolean pause = true;

// Font for text
PFont Font1 = createFont("Arial Bold", 18);

void setup() {
  // Resolution
  size (600, 600);

  // Instantiate matrices 
  A   = new float[width/cellSize][height/cellSize];
  A1  = new float[width/cellSize][height/cellSize];
  A2  = new float[width/cellSize][height/cellSize];
  A2n = new float[width/cellSize][height/cellSize];
  
  // Choose random initial parameters
  alpha = 1;
  beta  = random(0,1);
  gamma = random(0,.005);
alpha=1.00001;
beta=0.95;
gamma=0.001; 
  
  // Initialiaze max/min of A
  A_max = alpha;
  A_min = beta;

  // Disable antialiasing
  noSmooth();

  // Initialization of matrices
  for (int x=0; x<width/cellSize; x++) {
    for (int y=0; y<height/cellSize; y++) {
      A[x][y]   = beta;
      A2[x][y]  = beta;
      A2n[x][y] = beta;
      A1[x][y]  = 0;
    }
  }
  
  // Place seed in center
  A[(width/cellSize)/2][(width/cellSize)/2] = alpha;
  
  // Fill in background black
  background(0); 
}


void draw() {
  // Draw matrix
  for (int x=0; x<width/cellSize; x++) {
    for (int y=0; y<height/cellSize; y++) {
      float norm = (A[x][y]-A_min) * 1 / (A_max - A_min);
      fill(color(R * norm, G * norm, B * norm));
      rect (x*cellSize, y*cellSize, cellSize, cellSize);
    }
  }
  
  // Iterate if interval has been reached
  if (millis()-lastRecordedTime>interval) {
    if (!pause) {
      iteration();
      lastRecordedTime = millis();
    }
  }

  // Create new cells manually on pause
  if (pause && mousePressed) {
    // Map and avoid out of bound errors
    int xCellOver = int(map(mouseX, 0, width, 0, width/cellSize));
    xCellOver = constrain(xCellOver, 0, width/cellSize-1);
    int yCellOver = int(map(mouseY, 0, height, 0, height/cellSize));
    yCellOver = constrain(yCellOver, 0, height/cellSize-1);

    // Check against matrix
    if (A[xCellOver][yCellOver] == beta) {
      A[xCellOver][yCellOver] = alpha; // Place seed
    }
  } 
  
  // Output parameters and status
  textFont(Font1);
  fill(0, 255, 0);
  text("alpha = " + alpha, 10, 50);
  text("beta  = " + beta, 10, 70);
  text("gamma = " + gamma, 10, 90);
  if (!pause) {
    text("Simulating...", 10, 25);
  }
}

void iteration() {
  // Reset max/min values of A matrix
  A_max = Float.NEGATIVE_INFINITY;
  A_min = Float.POSITIVE_INFINITY;
  
  // Scan cells
  for (int x=1; x<width/cellSize - 1; x++) {
    for (int y=1; y<height/cellSize - 1; y++) {
      // Find min and max of the A matrix
      if (A[x][y] > A_max)
        A_max = A[x][y];
      else if (A[x][y] < A_min)
        A_min = A[x][y];
      if (x % 2 == 0) // odd row
      {
        // Grow ice
        if (A[x][y]>=alpha||A[x-1][y-1]>=alpha||A[x-1][y]>=alpha||A[x][y+1]>=alpha||
              A[x+1][y]>=alpha||A[x+1][y-1]>=alpha||A[x][y-1]>=alpha)
        {
          A1[x][y] = A[x][y] + gamma;
          A2[x][y] = 0;
        }
        else
        {
          A1[x][y] = 0;
          A2[x][y] = A[x][y];
        }
      }
      else // even row
      {
        // Grow ice
        if (A[x][y]>=alpha||A[x-1][y+1]>=alpha||A[x][y+1]>=alpha||A[x+1][y+1]>=alpha||
              A[x+1][y]>=alpha||A[x][y-1]>=alpha||A[x-1][y]>=alpha)
        {
          A1[x][y] = A[x][y] + gamma;
          A2[x][y] = 0;
        }
        else
        {
          A1[x][y] = 0;
          A2[x][y] = A[x][y];
        }
      }
      
    } // End of y loop
  } // End of x loop
  
  float avg_neigh;
  // Diffuse water
  for (int x=1; x<width/cellSize - 1; x++) {
    for (int y=1; y<height/cellSize - 1; y++) {
      if (x % 2 == 0)
        avg_neigh = (A2[x-1][y-1] + A2[x-1][y] + A2[x][y+1] + A2[x+1][y] + A2[x+1][y-1] + A2[x][y-1]) / 6;
      else
        avg_neigh = (A2[x-1][y+1] + A2[x][y+1] + A2[x+1][y+1] + A2[x+1][y] + A2[x][y-1] + A2[x-1][y]) / 6;
      A2n[x][y] = (A2[x][y] + avg_neigh) / 2;
    } // End of y loop
  } // End of x loop
  // Add updated water and ice
  for (int x=0; x<width/cellSize; x++) {
    for (int y=0; y<height/cellSize; y++) {
      A[x][y]  = A1[x][y] + A2n[x][y];
      A2[x][y] = A2n[x][y];
    }
  }
  //Print A
  for (int x=0; x<width/cellSize; x++) {
    for (int y=0; y<height/cellSize; y++) {
      print(A[x][y] + " ");
    }
    print("\n");
  }
  print("\n");
} // End of function

void keyPressed() {
  // R pressed - reset simulation
  if (key=='r') {
    setup();
    pause = true;
  }
  // SPACE BAR pressed - pause simulation
  if (key==' ') {
    pause = !pause;
  }
}

