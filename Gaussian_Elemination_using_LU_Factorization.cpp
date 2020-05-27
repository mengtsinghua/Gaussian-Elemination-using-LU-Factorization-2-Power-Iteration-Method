#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;

int main(int argc, char *argv[])
{ // FILE READINGS AND TAKING DATA FROM FILES
  string yaz; 
  int n = 0;  // n corresponds to row number of the matrix
  fstream A, b;
  A.open(argv[1]);//First parameter is the name of the file containing A matrix
  // A file is opened to get # of rows
  if (A.is_open())
  {
    while (getline(A,yaz))
    {
      n++; 
    }
    A.close();
  }else{
  	cout << "File couldn't open";
  }
  // 2-D matrix is created(A matrix)
  float **matrix = new float *[n];
  for (int i = 0; i < n; i++)
  {
    matrix[i] = new float[n];
  }
  // 1-D vector is created(b vector)
  float *vector = new float[n];
  // A.txt is opened to read the given A matrix and write into our 2-D created matrix
  string line;
  int sat1 = 0;
  A.open(argv[1]);
  while (getline(A, line))
  {
    istringstream reader(line);
    float dummy;
    for (int i = 0; i < n; i++)
    {
      reader >> dummy; // reading row by row
      matrix[sat1][i] = dummy; // writing element by element
    }
    sat1++;
  }
  // b.txt is opened to read the given b vector and write into our 1-D created vector
  b.open(argv[2]);//Second parameter is the name of the file containing b vector
  int sat2 = 0;
  while (getline(b, line))
  {
    istringstream reader(line);
    float dummy2;
    reader >> dummy2;
    vector[sat2] = dummy2;
    sat2++;
  }

  A.close();
  b.close();
  // checking if # of rows are matching each other. If not, then invalid input warning pops out 
  if (sat1 != sat2)
  {
    cout << "ERROR_1:INVALID INPUT.Matrix dimensions of A and b did not match!!! ";
  }
  else if (sat1 == sat2)
  {
    // Condition Number :if A matrix has 2x2 dimensions 
    if (n == 2) // (sat1==sat2==n just a reminder) 
    {
      int i, j = 0;
      float norm_1, norm_inf = 0; // norms are defined
      float sum_1, sum_inf = 0; // temporary variables are defined
      float det = (matrix[0][0] *matrix[1][1]) - (matrix[0][1] *matrix[1][0]);
      for (i = 0; i < n; i++)
      {
        for (j = 0; j < n; j++)
        {
          sum_inf = sum_inf + fabs(matrix[i][j]); // infinity norm adds absolute values of elements rowwise
          sum_1 = sum_1 + fabs(matrix[j][i]); // 1-norm adds absolute values of elements columnwise
        }
        if (sum_1 > norm_1)
        {
          norm_1 = sum_1; // if current column sum is greater than previous one, then change cond_1
        }
        if (sum_inf > norm_inf)
        {
          norm_inf = sum_inf; // if current row sum is greater than previous one, then change cond_inf
        }
        sum_1 = 0; // these variables are nullified  for every row or column cycle
        sum_inf = 0;
      }
      // since its a 2x2 matrix,no need to calculate inverse matrix.There occurs this relationship :
      float inv_norm_1 = norm_inf / fabs(det);// norms just swap but with 1/det coefficient
      float inv_norm_inf = norm_1 / fabs(det);
      // ConditionNumber=  norm at something x inverse norm at something(1,2,inf etc.)
      float cond_1 = norm_1 * inv_norm_1;
      float cond_inf = norm_inf * inv_norm_inf;

      cout << "The condition number of A matrix at norm 1 is : " << cond_1 << endl;
      cout << "The condition number of A matrix at norm infinity is : " << cond_inf << endl;
    }
    // ROW OPERATIONS
    
    int i, j, k, m = 0;
    int counter = 0;
    float temp[n] = {}; // temporary vector is created to use in row exchanges for A matrix
    float temp_b; // temporary varible is created to use in row exchanges for b vector
    
    // 1)Finding pivots and making row exchanges before any math operation

    for (j = 0; j < n; j++) // scanning column by column
    {
      for (i = 0; i < n; i++)
      {
        if (i == j) // After hitting a diagonal element :
        {
          float pivot = fabs(matrix[j][j]);// temporary pivot is created
          for (m = i + 1; m < n; m++)
          {//comparing the magnetidu between current pivot and elements under pivot's location in same column
            if (pivot < fabs(matrix[m][j]))
            {
              counter = m;// if a bigger pivot number exists, then register its row location
              pivot = fabs(matrix[m][j]);// change previous pivot with a bigger current one to cycle through
            }
          }
          //After getting row location of possible biggest pivot in same column under the diagonal:
          for (k = 0; k < n; k++)
          {
            if (counter != 0)// Make row exchanges between i'th row and counter'th row if necessary
            { 
              temp[k] = matrix[counter][k]; // keep row elements of counter'th row
              matrix[counter][k] = matrix[i][k];// swap down the i'th row to the counter'th row
              matrix[i][k] = temp[k]; // the row containing biggest abs(number) in same column under the diagonal
			  //becomes the row containing pivot which will have location in the diagonal  
              if (k == 0)
              {//Apply same swapping rule for b vector except just for once therefore k==0
                temp_b = vector[counter];
                vector[counter] = vector[i];
                vector[i] = temp_b;
              }
            }
          }
        } // i==j
        counter = 0; // Counter is refreshed for every cycle
      } // i loop
    } // j loop
    // 2)After alligning the pivots in the diagonal, its time to apply gaussian elimination
    int row, column = 0;
    float multiplier;
    for (int i = 0; i < n; i++)//this loop exists to prevent applying elimination in pivot's row 
    {
      for (int row = 0; row < n; row++) //scanning row by row 
      {
        if (row != i)//If the row is not pivot's row, then apply gaussian elimination
        {// the value of element which has same column index with current pivot/the value of pivot equals to a temporary multiplier
          multiplier = matrix[row][i] / matrix[i][i]; 
          //Row operation for b vector but just for once since its 1-D
          vector[row] = vector[row] - (multiplier * vector[i]);
          for (int column = 0; column < n; column++)
          {// substract each element in a non-pivot row with this multiplier times corresponding element in the pivot row
            matrix[row][column] = matrix[row][column] - (multiplier * matrix[i][column]);
          }
        }
      }
    }
    // Machine Precision Check
    // Make both matrix and vector elements zero if the abs value of element smaller than epsilon
    for (i = 0; i < n; i++)
    {
      if (fabs(vector[i]) < 3.49209e-006)
      { // it was 1,1920929e-007 for IEEE SP but i narrowed down the precision by choice
        vector[i] = 0;
      }
      for (j = 0; j < n; j++)
      {
        if (fabs(matrix[i][j]) < 3.49209e-006)
        {
          matrix[i][j] = 0;
        }
      }
    }
    // Singularity Check
    int singular = 0;// This variable added to print out ERROR_2 just for once in case there may have multiple zero rows
    int count = 0;
    
    for (i = 0; i < n; i++)// scanning row by row
    {
      count = 0;
      for (j = 0; j < n; j++)
      {//If an element in a row equals to zero or nan, then increase count by one
        if (matrix[i][j] == 0 || (matrix[i][j] != matrix[i][j]))
        {
          count++;
        }
        if (count == n && singular == 0)// If all elements are zero, then the matrix is singular 
        {
          cout << "ERROR_2:The matrix is singular!!!" << endl;
          singular++;
          break;
        }
      }
    }
    // Creating Solution Vector:x
    float *solution;
    solution = new float[n];
    // A matrix is now a diagonal matrix,therefore Ax=b relationship can easily be implented :
    for (i = 0; i < n; i++)
    {// Back substitution
      solution[i] = (vector[i] / matrix[i][i]);
    }
    // Writing into x.txt file
    if (singular == 0)//If matrix is not singular, then write into x.txt.If not, then there is no need to write into x.txt
    {
      ofstream sol_x("x.txt");
      if (sol_x.is_open())
      {
        for (int i = 0; i < n; i++)
        {
          sol_x << solution[i] << endl;//Transfer the solution matrix into x.txt file
        }
        sol_x.close();
      }
    }
  }
  //Delete dynamic arrays to free memory
  delete [] matrix;
  delete [] vector;
  
  return 0;
}
