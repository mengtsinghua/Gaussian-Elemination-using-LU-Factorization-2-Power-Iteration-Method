#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;


class Matrix {
	public:
	int r,c;
	float **array;
	float *vec;
	Matrix(){};
	Matrix(int,int);
	Matrix(int,int,float *);
	Matrix(int,int,float **);
	Matrix operator * (float);
	Matrix operator / (float);
	
};

Matrix::Matrix(int row,int column){
     r=row;
     c=column;
     array = new float *[r];
     vec= new float [r];
    for(int i = 0; i < r; i++){
    	  array[i] = new float [c];
    	  vec[i]=0;
	}
	for(int i=0;i<r;i++){
		for(int j=0;j<c;j++){
			array[i][j]=0;
		}
	}
	cout<<"error1"<<endl;
} 

Matrix::Matrix(int row,int column,float ** pass){
	r=row;
	c=column; 
    array=pass;
} 

Matrix::Matrix(int row,int column,float * pass){
	r=row;
	c=column; 
    vec=pass;
} 




Matrix product_1(Matrix first,Matrix second){ // MATRIX*MATRIX
	Matrix temp(first.r,second.c);
     for(int i = 0; i < first.r; i++){
        for(int j = 0; j < second.c; j++){
            for(int k = 0; k < first.c; k++)
            {
             temp.array[i][j] += first.array[i][k] * second.array[k][j];
            }}}
           cout<<"error2"<<endl;
    return temp;
    
}
Matrix product_2(Matrix first,Matrix second){ //MATRIX*VECTOR
	Matrix temp(first.r,1);
	for(int i = 0; i < first.r; i++){
        for(int j = 0; j < first.c; j++){
           
             temp.vec[i] += first.array[i][j] * second.vec[j];
            }}
            cout<<"error3"<<endl;
	return temp;
}
Matrix product_3(Matrix first,Matrix second){ // VECTOR*VECTOR=MATRIX

	Matrix temp(first.r,second.c);
        // 4*1 by 1*4 
		for(int i = 0; i < first.r; i++){
        for(int j = 0; j < second.c; j++){
           
             temp.array[i][j] += first.vec[i] * second.vec[j];
            }}
            cout<<"error4"<<endl;
            return temp;
	
	
}
float number(Matrix first,Matrix second){ // VECTOR*VECTOR=FLOAT
		float result=0;
	 // 1*4 by 4*1
		for(int i = 0; i < second.r; i++){
            result += first.vec[i] * second.vec[i];
            }
            cout<<"error5"<<endl;
            return result;
	
}

float norm_finder (Matrix v){
  float	max=0;
	for(int i=0;i<v.r;i++){
		if(fabs(v.vec[i])>max){
			max=fabs(v.vec[i]);
		}	
	}
	cout<<"error6"<<endl;
return max;
}

float magnitude(Matrix v){
	float sum=0;
	for(int i=0;i<v.r;i++){
		sum += (v.vec[i]*v.vec[i]);
	}
	cout<<"error7"<<endl;
	return sqrt(sum);
}

Matrix transpose_1(Matrix H){ // MATRIX TRANSPOSE
	Matrix temp(H.c,H.r);
	for(int i=0;i<H.r;i++){
		for(int j=0;j<H.c;j++){
			temp.array[j][i]=H.array[i][j];
		}
	}
	cout<<"error8"<<endl;
	return temp;
}
Matrix transpose_2(Matrix H){ // VECTOR TRANSPOSE
	Matrix temp(H.c,H.r);
	for(int i=0;i<H.r;i++){
		
			temp.vec[i]=H.vec[i];
		
	}
	cout<<"error9"<<endl;
	return temp;
}

Matrix Matrix::operator*(float multiplier) {  
       for (int i=0;i<r;i++) {
        vec[i] = vec[i]*multiplier;
		}
		cout<<"erro10"<<endl;
       return *this;
	   }
Matrix Matrix::operator / (float divider) {  
       for (int i=0;i<r;i++) {
        vec[i] = vec[i]/divider;
		}
		cout<<"error11"<<endl;
       return *this;
	   }	   
	   












int main()
{

  string yaz; 
  int n = 0;  // n corresponds to row number of the matrix
  fstream A, b;
  A.open("A.txt");//First parameter is the name of the file containing A matrix
  // A file is opened to get # of rows
  if (A.is_open())
  {
    while (getline(A,yaz))
    {
      n++; 
    }
    cout<<"error12"<<endl;
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

 
  // A.txt is opened to read the given A matrix and write into our 2-D created matrix
  string line;
  int sat1 = 0;
  A.open("A.txt");
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

/*
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
	cout << matrix[i][j] << " ";
	}
	cout << endl;
}
*/
  A.close();

float *vector = new float[n];
for(int i=0;i<n;i++){
	vector[i]=1;
}

Matrix data(n,n,matrix);
Matrix initial_vec(n,1,vector);
Matrix ret(n,1);
//ret=product_2(data,initial_vec);
//product_2(data,initial_vec);
//cout << ret.vec[1];
//cout << data.array[0][0];
//cout << norm_finder(initial_vec);
//cout << magnitude(initial_vec);
//cout << transpose(data).array[0][1];




float current_norm=0;
float prev_norm=0;
float tolerance=100;
Matrix temp(n,1);
cout<<"error13"<<endl;
while(tolerance>1e-6){
	
	temp=product_2(data,initial_vec);
	prev_norm=current_norm;
	current_norm=norm_finder(temp);
	temp=temp/current_norm;
	initial_vec=temp;
    tolerance=fabs(current_norm-prev_norm);
	
}
cout<<"error14"<<endl;

/*
cout << current_norm << endl;
float tests=0;
for(int i=0;i<initial_vec.r;i++){

	tests += initial_vec.vec[i]*initial_vec.vec[i];
	
	cout << tests << endl;
}
*/

// Householder - Deflation

Matrix first_element(n,1);
first_element.vec[0]=1;
transpose_2(initial_vec);
//cout << product_3(initial_vec,transpose_2(initial_vec)).array[3][3];

//cout << product_2(transpose_1(initial_vec),first_element);




//float bek=number(transpose_2(initial_vec),initial_vec);
//cout << bek ;
Matrix v(n,1);
for(int i=0;i<n;i++){
      v.vec[i]	=initial_vec.vec[i]-(first_element*magnitude(initial_vec)).vec[i];
    // cout << v.vec[i] << endl;
	 }
 // eye matrix
 Matrix eye(n,n);
 for(int i=0;i<n;i++){
 	for(int j=0;j<n;j++){
 		if(i==j){
 			eye.array[i][j]=1;
		 }else{
		 	eye.array[i][j]=0;
		 }
	 }
 }
cout<<"error15"<<endl;

Matrix temp_m(n,n);
temp_m=product_3(v,transpose_2(v));
//temp_m=temp_m*2;
float temp_n=number(transpose_2(v),v);

for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		float temp = temp_m.array[i][j]*2;
		temp_m.array[i][j]=temp/temp_n;
	}
}

Matrix Householder(n,n);

for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		Householder.array[i][j]=eye.array[i][j]-temp_m.array[i][j];
		
	}
	
}

cout<<"error16"<<endl;
Matrix Householder_trans(n,n);

Householder_trans=transpose_1(Householder);


Matrix Z(n,n);

Z = product_1(product_1(Householder,data),Householder_trans);




 //  Matrix temp_s(n,n);
  // temp_s= (temp_m /temp_n);
   

for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
	if(fabs(Z.array[i][j])<2e-6){
		Z.array[i][j]=0;
	}else{
		continue;
	}
	}
	
}

Matrix Z_trim(n-1,n-1);

for(int i=0;i<n-1;i++){
	for(int j=0;j<n-1;j++){
	
		Z_trim.array[i][j]=Z.array[i+1][j+1];
		};
	}


float *vector_2 = new float[n-1];
for(int i=0;i<n-1;i++){
	vector_2[i]=1;
}


cout<<"error17"<<endl;
Matrix initial_vec_2(n-1,1,vector_2);


float current_norm_2=0;
float prev_norm_2=0;
float tolerance_2=100;
Matrix temp_2(n-1,1);

while(tolerance_2>120){
	
	
	temp_2=product_2(Z_trim,initial_vec_2);
	prev_norm_2=current_norm_2;
	current_norm_2=norm_finder(temp_2);
	temp_2=temp_2/current_norm_2;
	initial_vec_2=temp_2;
    tolerance_2=fabs(current_norm_2-prev_norm_2);
	
}


cout<<"error18"<<endl;
cout << current_norm << " "<< initial_vec.vec[2];
/*
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
	
	cout << 	Z.array[i][j] <<" ";
	
	}
	cout << endl;
}
*/


return 0;
}
