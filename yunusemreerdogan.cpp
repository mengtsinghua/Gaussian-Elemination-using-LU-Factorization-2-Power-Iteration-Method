#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;

class Matrix
{
  public:
    int r, c; // row and column
    float **array; // Matrix or vector construction inside the class

    Matrix(){}; // default constructor
    Matrix(int, int); // constructor for a new empty matrix or vector
    Matrix(int, int, float **); // constructor for a known given matrix
    Matrix operator *(float); // matrix multiplication with a scalar number 
    Matrix operator / (float); // matrix division with a scalar number

};
Matrix::Matrix(int row, int column)
{
  r = row;
  c = column;
  array = new float *[r];
  for (int i = 0; i < r; i++)
  {
    array[i] = new float[c];
  }
}


Matrix::Matrix(int row, int column, float **pass)
{
  // passing A matrix into class constructor
  r = row;
  c = column;
  array = new float *[r];
  for (int i = 0; i < r; i++)
  {
    array[i] = new float[c];
  }
  array = pass;
}

Matrix product_1(Matrix first, Matrix second)
{
  // MATRIX*MATRIX multiplication
  Matrix temp(first.r, second.c); 
  // temporary object is created to store the result as an object
  for (int i = 0; i < first.r; i++)
  {
    for (int j = 0; j < second.c; j++)
    {
      for (int k = 0; k < first.c; k++)
      {
        temp.array[i][j] += first.array[i][k] *second.array[k][j];
      }
    }
  }

  return temp;
}

float norm_finder(Matrix v)
{
  // Finds the largest absolute value in a vector column
  float max = 0;
  for (int i = 0; i < v.r; i++)
  {

    if (fabs(v.array[i][0]) > max)
    {
      max = fabs(v.array[i][0]);
    }
  }
  return max;
}

float magnitude(Matrix v)
{
  // 2-norm of a vector ,i.e. geometric magnitude c^2=a^2+b^2;
  float sum = 0;
  for (int i = 0; i < v.r; i++)
  {
    sum += (v.array[i][0] *v.array[i][0]);
  }
  return sqrt(sum);
}

Matrix transpose_1(Matrix H)
{
  // MATRIX TRANSPOSE
  Matrix temp(H.c, H.r);
  for (int i = 0; i < H.r; i++)
  {
    for (int j = 0; j < H.c; j++)
    {
      temp.array[j][i] = H.array[i][j];
    }
  }
  return temp;
}


Matrix Matrix::operator *(float multiplier)
{
  // Matrix * scaler number
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < c; j++)
    {
      array[i][j] = array[i][j] *multiplier;
    }
  }
  return  *this;
}

Matrix Matrix::operator / (float divider)
{
  // Matrix / scaler number
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < c; j++)
    {
      array[i][j] = array[i][j] / divider;
    }
  }
  return  *this;
}

int main(int argc, char *argv[])
{
  // Reading from A.txt file
  string yaz;
  int n = 0; // n corresponds to row number of the matrix
  fstream A;
  A.open(argv[1]); //First parameter is the name of the file containing A matrix
  // A file is opened to get # of rows
  if (A.is_open())
  {
    while (getline(A, yaz))
    {
      n++;
    }

    A.close();
  }
  else
  {
    cout << "File couldn't open";
  }
  // 2-D matrix is created(A matrix)
  float **matrix = new float *[n];
  for (int i = 0; i < n; i++)
  {
    matrix[i] = new float[n];
  }
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
  A.close();


  // Write given matrix into data object
  Matrix data(n, n, matrix);
  // Create an empty object with 1-D vector
  Matrix initial_vec(n, 1);

  // Initialized vector = [1,1,1,......1] ;
  for (int i = 0; i < initial_vec.r; i++)
  {
    initial_vec.array[i][0] = 1;
  }

  float current_norm = 0; //starting eigenvalue
  float prev_norm = 0; //previous eigenvalue
  float tolerance = 100; //a convenient value to start the while loop
  Matrix temp(n, 1); // temporary vector object created
  float prev_tolerance = 0; // previous difference between current_norm and prev_norm
  int sign_control = 0; //checking whether the sign of the first element of initial_vec is changed or not 
  int counter = 0; // counts the iteration 

  while (tolerance > atof(argv[2]))
  {
    counter++;
    temp = product_1(data, initial_vec); // A * x(initial vector)
    prev_norm = current_norm; // updating eigenvalue
    current_norm = norm_finder(temp); // normalize the eigenvalue 

    prev_tolerance = tolerance; // updating the previous tolerance

    temp = temp / current_norm; // normalize the x vector

    if (temp.array[0][0] *initial_vec.array[0][0] < 0)
    {
      sign_control++;// if the sign of the first element of initial_vec is changed,then add 1 
      // this process exists to detect the negative largest eigenvalue
    }

    initial_vec = temp; // updating the x vector
    tolerance = fabs(current_norm - prev_norm); // updating the current tolerance
    
	if (prev_tolerance - tolerance == 0)
    {
      break; // if the latest tolerance difference is equal, no need to restart the loop
    }
  }
  
  // if the sign change occurs more than half of the number of iterations,then change sign of the eigenvalue
  if (sign_control >= counter / 2)
  {
    current_norm =  - current_norm;
  }
  // Householder - Deflation

  // To obtain second largest eigenvalue, we need to omit the first row and column of H*A*H' matrix
  // H represents Householder matrix

  //the equations : v=initial_vec-(norm(initial_vec,2))*first_element;
  //the equations : Householder(H)=eye(n)-(2*v*v')/(v'*v); eye is identity matrix with n dimensions
  //the equations : Z=H*A*H'
  //the equations : Z_trim = Z(second row : n ; second column : n ); first row and column are omitted

  Matrix first_element(n, 1);
  first_element.array[0][0] = 1; // only the first element is 1, others 0

  Matrix v(n, 1); // Create v vector
  for (int i = 0; i < n; i++)
  {
    v.array[i][0] = initial_vec.array[i][0] - (first_element *magnitude(initial_vec)).array[i][0];
  }
  // eye matrix i.e. identity matrix
  Matrix eye(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j)
      {
        eye.array[i][j] = 1;
      }
      else
      {
        eye.array[i][j] = 0;
      }
    }
  }

  Matrix temp_m(n, n); // temporay object to store v*v' ;
  temp_m = product_1(v, transpose_1(v));

  Matrix number(1, 1); // temporary object to store v'*v ;
  number = product_1(transpose_1(v), v);

  // calculation of (2*v*v')/(v'*v)
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      float temp = temp_m.array[i][j] *2;
      temp_m.array[i][j] = temp / number.array[0][0];
    }
  }

  Matrix Householder(n, n); // Create H matrix

  // calculation of Householder(H)=eye(n)-(2*v*v')/(v'*v);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      Householder.array[i][j] = eye.array[i][j] - temp_m.array[i][j];
    }
  }

  Matrix Householder_trans(n, n);
  Householder_trans = transpose_1(Householder); // H transpose is created
  Matrix Z(n, n);
  Z = product_1(product_1(Householder, data), Householder_trans); // Z=H*A*H' 

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (fabs(Z.array[i][j]) < 2e-6)
      {
        Z.array[i][j] = 0;
      }
      else
      {
        continue;
      }
    }
  }

  Matrix Z_trim(n - 1, n - 1); // trim matrix is created

  for (int i = 0; i < n - 1; i++)
  {
    for (int j = 0; j < n - 1; j++)
    {
      Z_trim.array[i][j] = Z.array[i + 1][j + 1]; // one dimension of Z matrix is omitted
    };
  }

  Matrix initial_vec_2(n - 1, 1); // new x vector is created 
  for (int i = 0; i < initial_vec_2.r; i++)
  {
    initial_vec_2.array[i][0] = 1; // vector initialized with ones
  }


  // EXACTLY SAME PROCESS WITH PREVIOUS ITERATION LOOP

  float current_norm_2 = 0;
  float prev_norm_2 = 0;
  float tolerance_2 = 100;
  Matrix temp_2(n - 1, 1);
  float prev_tolerance_2 = 0;

  int sign_control_2 = 0;
  int counter_2 = 0;

  while (tolerance_2 > atof(argv[2]))
  {
    counter_2++;
    temp_2 = product_1(Z_trim, initial_vec_2);
    prev_norm_2 = current_norm_2;
    current_norm_2 = norm_finder(temp_2);

    prev_tolerance_2 = tolerance_2;

    temp_2 = temp_2 / current_norm_2;

    if (temp_2.array[0][0] *initial_vec_2.array[0][0] < 0)
    {
      sign_control_2++;
    }

    initial_vec_2 = temp_2;
    tolerance_2 = fabs(current_norm_2 - prev_norm_2);

    if (prev_tolerance_2 - tolerance_2 == 0)
    {
      break;
    }

  }

  if (sign_control_2 >= counter_2 / 2)
  { //Sign of second eigenvalue is also checked
    current_norm_2 =  - current_norm_2;
  }

  // IF NOT ENOUGH ARGUMENTS THEN ERROR OCCURS
  if (argc != 4)
  {
    cout << "ERROR DETECTED:INVALID # OF ARGUMENTS.Program Could Not Initialized.";
  }
  else
  {
    // IF NO ERRORS, THEN PRINT AND WRITE THE OUTPUT
    cout << "Eigenvalue#1: " << current_norm << endl;

    for (int i = 0; i < n; i++)
    {
      cout << initial_vec.array[i][0] << " ";
      cout << endl;
    }
    cout << "Eigenvalue#2: " << current_norm_2 << endl;
    ofstream sol_x(argv[3]);
    if (sol_x.is_open())
    {
      sol_x << "Eigenvalue#1: " << current_norm << endl;
      for (int i = 0; i < n; i++)
      {
        sol_x << initial_vec.array[i][0] << endl; 
        //Transfer the eigenvector into x.txt file
      }
      sol_x << "Eigenvalue#2: " << current_norm_2 << endl;
      sol_x.close();
    }
  }

  // DELETE MAINLY USED DYNAMIC MEMORY ARRAYS
  delete [] data.array;
  delete [] Z_trim.array;
  delete [] initial_vec.array;
  delete [] initial_vec_2.array;

  return 0;
}
