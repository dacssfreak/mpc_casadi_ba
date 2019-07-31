/*
 * Copyright (C) 2013- Acho Arnold
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "matrix.h"
#include <windows.h>

int ALLOC;

void matrix_print (Matrix *matrix)
{
	int i,j;
	printf("\n");
	for(i = 0; i < matrix->row_size;i++)
	{
		printf("\t\t");
		for(j = 0; j < matrix->col_size;j++)
		{
			printf("%9.2f", matrix->matrix_entry[i][j]);
		}
   	printf("\n");
	}
	printf("\n");
}

void matrix_print_detail(Matrix* matrix, char* name, int rep) {
	printf("\n %s_%d is:\n", name, rep);
	matrix_print(matrix);
}

void matrix_print_part(Matrix *matrix, int start_index)
{
    int j,i;
    for (i = 0; i < matrix->row_size; ++i)\
    {
        for (j = start_index; j < matrix->col_size; ++j)
        {
            printf("\t\t%9.2f", matrix->matrix_entry[i][j]);
        }
	printf("\n");
    }
}

void matrix_fill(Matrix *matrix)
{
	int i,j;
	printf("Enter the contents of the matrix (%d/%d):\n", matrix->row_size, matrix->col_size);
	for(i=0; i < matrix->row_size; i++)
	{
		for(j=0;j < matrix->col_size; j++)
		{
			scanf_s("%f",&matrix->matrix_entry[i][j]);
		}
	}
}

/*Function to create an identity matrix	*/
Matrix * matrix_callalloc(int matrix_size)
{
	Matrix *result = matrix_alloc(matrix_size, matrix_size);
	int i,j;
	
	for (i = 0; i < matrix_size; i += 1)
	{
		for (j = 0; j < matrix_size; j += 1)
		{
			if (j == i)
			{
				result->matrix_entry[i][j] = 1;
			}
			
			else
			{
				result->matrix_entry[i][j] = 0;
			}
		}
	}
	
	return result;
}

Matrix *matrix_alloc(int row_size, int col_size)
{
	int j;
	Matrix* new_matrix = (Matrix*)malloc(sizeof(Matrix));

	if (new_matrix != NULL) {
		//Allocating memory for the new matrix structure
		new_matrix->row_size = row_size;
		new_matrix->col_size = col_size;
		size_t size = new_matrix->row_size * sizeof(float*);
		new_matrix->matrix_entry = malloc(size);
		if (new_matrix->matrix_entry) {
			for (j = 0; j < new_matrix->row_size; j++)
			{
				float f;
				size_t size2 = new_matrix->col_size * sizeof(f);
				new_matrix->matrix_entry[j] = malloc(size2);
			}
		}
	}
	ALLOC = ALLOC + 1;
	return new_matrix;
}

/*Copies Matrix1 into matrix2 */
void matrix_copy(Matrix *matrix1, Matrix *matrix2)
{
	int i, j;
	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (j = 0; j < matrix1->col_size; j += 1)
		{
			matrix2->matrix_entry[i][j] = matrix1->matrix_entry[i][j];
		}
	}
}


Matrix* matrix_multiply(const Matrix *matrix1, const Matrix *matrix2)
{
	int i, j, k;
	float sum;

	if (matrix1->col_size != matrix2->row_size)
	{
		terminate("ERROR: The number columns of matrix1  != number of rows in matrix2!");
	}

	Matrix* result = matrix_alloc( matrix1->row_size,matrix2->col_size);

	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (k = 0; k < matrix2->col_size; k += 1)
		{
			sum = 0;
			
			for (j = 0; j < matrix1->col_size; j += 1)
			{
				sum += matrix1->matrix_entry[i][j] * matrix2->matrix_entry[j][k];
			}
			
			result->matrix_entry[i][k] = sum;
		}
	}

    return result;
}

Matrix* matrix_pow(Matrix *matrix, int index)
{
	if(index == 1)
	{
	  Matrix  *result = matrix_alloc (matrix->row_size, matrix->col_size);
	  matrix_copy(matrix, result);
	  return result;
	}
	else
	{
	
		int i, j,k,l,count;
		float sum;
		
		Matrix  *temp = matrix_alloc (matrix->row_size, matrix->col_size); //Allocating space for a temporal matrix
		Matrix  *result = matrix_alloc (matrix->row_size, matrix->col_size); //Allocating space for the result matrix
	
		matrix_copy(matrix, temp);
		
		count = index/2 -1;
		if (count < 1)
		{
			matrix_copy(matrix, result);
		}
		
		else
		{
			for (l = 0; l < count; l += 1)
			{
				for (i = 0; i < matrix->row_size; i += 1)
				{
					for (k = 0; k < matrix->col_size; k += 1)
					{
						sum = 0;
			
						for (j = 0; j < matrix->col_size; j += 1)
						{
							sum += (temp->matrix_entry[i][j] * matrix->matrix_entry[j][k]);
						}
					
						result->matrix_entry[i][k] = sum;
					}
				}
			
				/* Copying the result matrix into the temp matrix for further 
				 * multiplication */
				matrix_copy(result, temp);
			}
		}

		/*	Freeing the temp matrix		*/
		matrix_free(temp);
		if (index%2 == 0)
		{
			Matrix *result_final = matrix_multiply(result, result);
			/* Freeing the result Matrix	*/
			matrix_free(result);
		
			return result_final;
		}
		
		else
		{
			Matrix *temp = matrix_multiply(matrix, result);
			Matrix *result_final = matrix_multiply(temp, result);
			
			/* Freeing the temp matrix		*/
			matrix_free(temp);
			
			/* Freeing the result Matrix	*/
			matrix_free(result);
		
			return result_final;
		}//End of else statement
	}
}

void matrix_free(Matrix *matrix)
{
	int j;
	for (j = 0; j < matrix->row_size; j++)
	{
		free(matrix->matrix_entry[j]);
	}
	free(matrix->matrix_entry);
	free(matrix);
	ALLOC = ALLOC - 1;
}

void print_alloc()
{
	printf("\n ALLOC is: %d", ALLOC);
}

/*Function which divides all row entries by the value of a the diagonal */
void row_divide(Matrix *matrix, int pivot)
{
    int j;
    float 	divisor = matrix->matrix_entry[pivot][pivot], 
              result;

    for(j = pivot; j < matrix->col_size; j++)
    {
              result = (matrix->matrix_entry[pivot][j]  /  divisor);
              matrix->matrix_entry[pivot][j] = result;
    }

}

 /*Function to carry out row operations*/
void row_operation(Matrix *multiplier_matrix,Matrix *matrix, int pivot, int row_index)
{
    int j;
    float multiplier = (matrix->matrix_entry[row_index][pivot] / matrix->matrix_entry[pivot][pivot]);
    //Loop which checks if matrix is provided to store the multiplier
    if(multiplier_matrix != NULL)
      {
	multiplier_matrix ->matrix_entry[row_index][pivot] = multiplier;
      }

    
    for(j=0; j < matrix->col_size; j++)
    {
	    matrix->matrix_entry[row_index][j] -=  multiplier * matrix->matrix_entry[pivot][j];
    }
}

void matrix_row_reduce( Matrix *matrix, int zero_control )
{
    int pivot, row_index;
    float multiplier;
    for( pivot = 0; pivot < matrix->row_size ; pivot++)
    {
         
      error_zeros(matrix, zero_control); //Function checks if there are too many zeros in a single row
	    if(	(matrix->matrix_entry[pivot][pivot] != 1) && (matrix->matrix_entry[pivot][pivot] != 0)	)
	    {
		row_divide(matrix, pivot);
	    }

	    for (row_index = pivot+1; row_index < matrix->row_size; row_index++)
	    {
		    if (matrix->matrix_entry[pivot][pivot] != 0)
		    {
		      row_operation(NULL,matrix, pivot, row_index);
		    }
	    }

		for(row_index = pivot-1; row_index >=0; row_index --)
		{
			if (matrix->matrix_entry[pivot][pivot] != 0)
			{
			  row_operation(NULL,matrix, pivot, row_index);
			}
		}
	}
}

void LU_decompose(Matrix *upper_triangular, Matrix *lower_triangular)
{
 int pivot, row_index;
    float multiplier;
    for( pivot = 0; pivot < upper_triangular->row_size ; pivot++)
    {
         
      error_zeros(upper_triangular, upper_triangular->col_size); //Function checks if there are too many zeros in a single row
	    for (row_index = pivot+1; row_index < upper_triangular->row_size; row_index++)
	    {
		    if ( upper_triangular->matrix_entry[pivot][pivot] != 0)
		    {

		      row_operation(lower_triangular,upper_triangular, pivot, row_index);
		    }
	    }
    }
}


Matrix* matrix_subtract(Matrix* matrix1, Matrix* matrix2)
{
	int i, j;

	if (!(matrix_equal_size(matrix1, matrix2)))
	{
		terminate("ERROR: The matrices you are trying to subtract have different sizes");
	}
	Matrix* result = matrix_alloc(matrix1->row_size, matrix1->col_size);
	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (j = 0; j < matrix1->col_size; j += 1)
		{
			result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] - matrix2->matrix_entry[i][j];
		}
	}

	return result;
}


Matrix* matrix_add(Matrix* matrix1, Matrix* matrix2)
{
	int i, j;

	if (!(matrix_equal_size(matrix1, matrix2)))
	{
		terminate("ERROR: The matrices you are trying to add have different sizes");
	}
	Matrix* result = matrix_alloc(matrix1->row_size, matrix1->col_size);
	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (j = 0; j < matrix1->col_size; j += 1)
		{
			result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] + matrix2->matrix_entry[i][j];
		}
	}

	return result;
}


int matrix_equal_size( Matrix *matrix1, Matrix *matrix2)
{

  return (matrix1->row_size == matrix2->row_size && \
		  matrix1->col_size == matrix2->col_size);
}



/*
  This function checks if there is a line containing too many zero's and it exits
  if such a line is found
*/
void error_zeros( Matrix *matrix, int control_index)
{
      int i,j,count;
      for(i=0; i<matrix->row_size; i++)
      {
	    count=0;
	    for(j = 0;  j < matrix->col_size; j++)
	    {
	      if( matrix->matrix_entry[i][j] == 0)
	      {
		count++;     
	      }
	      else
	      {
		return;
	      }
	      if(count == control_index)
	      {
		fprintf(stdout,"\nProcess fail because row %d contains %d  zeros\n",i+1,control_index);
		matrix_print(matrix);
		exit(1);
	      }
	    }
	  }
}  


void terminate (char * string)
{
  fprintf(stdout,"\n%s\n",string);
  fprintf(stdout,"The program is exiting now. . . .\n\n");
  exit(-1);
}

//  Inverse einer Matrix
//  Berechung nach dem Gauss-Jordan Verfahren
//  Lit. Helmut Selder, Einführung in die numerische Mathematik für Ingenieure, HANSER
//
//  Das Verfahren löst die n Gleichungssysteme (für je eine Spalte der Einheitsvektoren) 
//  in einem gemeinsamen Eliminationsverfahren.
//  Statt nur einem Vektor für die rechte Seite, erweitert man die Matrix um alle n Einheitsvektoren
//  a11 a12 a13  | 1  0  0
//  a21 a22 a23  | 0  1  0
//  a31 a32 a33  | 0  0  1
//  und eliminiert in der Matrix a alle Elemente ausserhalb der Diagonalen
//  zusätzlich sorgt man durch eine Division fÜr lauter 1 in der Diagonalen
//  Das schaut dann so aus:
//   x1  x2  x2
//   1   0   0   | b11 b12 b13
//   0   1   0   | b21 b22 b23
//   0   0   1   | b31 b32 b33
//  Dadurch fällt die Berechnung der Lösungen aus der Dreiecksmatrix (xn, xn-1, ... x1)
//  weg, weil die Lösungen sofort abgelesen werden können.
//  x1 = b11  bzw. b12 oder b13
//  x2 = b21  bzw. b22 oder b23
//  x3 = b31  bzw. b23 oder b33
//  Das bedeutet, die Einheitsmatrix ist bei der Elimination in die Matrix B übergegangen,
//  und das ist die Inverse Matrix.

//Hint: a will be modified
int matrix_inverse(Matrix* a, Matrix* ainv)
{
	int n = a -> col_size;
	int   i, j;                    // Zeile, Spalte
	int   s;                       // Elimininationsschritt
	int   pzeile;                  // Pivotzeile
	int   fehler = 0;              // Fehlerflag
	float f;                      // Multiplikationsfaktor
	const double Epsilon = 0.01;   // Genauigkeit
	float Maximum;                // Zeilenpivotisierung
	int pivot = 1;

	for (int k = 0; k < a->row_size; k++) {

		float* as = a->matrix_entry[k];
		a->matrix_entry[k] = malloc(a->col_size * sizeof(float) * 2);

		if (a->matrix_entry[k])
			for (int x = 0; x < a->col_size; x++)
				a->matrix_entry[k][x] = as[x];
		free(as);
	}
	a->col_size *= 2;

	// ergänze die Matrix a um eine Einheitsmatrix (rechts anhängen)
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
		{
			a->matrix_entry[i][n + j] = 0.0;
			if (i == j)
				a->matrix_entry[i][n + j] = 1.0;
		}
	}

	// die einzelnen Eliminationsschritte
	s = 0;
	do {
		// Pivotisierung vermeidet unnötigen Abbruch bei einer Null in der Diagnonalen und
		  // erhöht die Rechengenauigkeit
		
		Maximum = myNewAbsAreSoHard(a->matrix_entry[s][s]);
		
		if (pivot)
		{
			pzeile = s;
			for (i = s + 1; i < n; i++)
				if (myNewAbsAreSoHard((double)a->matrix_entry[i][s]) > Maximum)
				{
					Maximum = myNewAbsAreSoHard((double)a->matrix_entry[i][s]);
					pzeile = i;
				}
		}
		fehler = (Maximum < Epsilon);

		if (fehler) {
			fehler = 1;
			break;           // nicht lösbar 
		}


		if (pivot)
		{
			if (pzeile != s)  // falls erforderlich, Zeilen tauschen
			{
				float h;
				for (j = s; j < 2 * n; j++) {
					h = a->matrix_entry[s][j];
					a->matrix_entry[s][j] = a->matrix_entry[pzeile][j];
					a->matrix_entry[pzeile][j] = h;
				}
			}
		}

		// Eliminationszeile durch Pivot-Koeffizienten f = a->matrix_entry[s][s] dividieren
		f = a->matrix_entry[s][s];
		for (j = s; j < 2 * n; j++)
			a->matrix_entry[s][j] = a->matrix_entry[s][j] / f;

		// Elimination --> Nullen in Spalte s oberhalb und unterhalb der Diagonalen
		// durch Addition der mit f multiplizierten Zeile s zur jeweiligen Zeile i
		for (i = 0; i < n; i++) {
			if (i != s)
			{
				f = -a->matrix_entry[i][s];                 // Multiplikationsfaktor
				for (j = s; j < 2 * n; j++)    // die einzelnen Spalten
					a->matrix_entry[i][j] += f * a->matrix_entry[s][j];       // Addition der Zeilen i, s
			}
		}

		s++;
	} while (s < n);

	if (fehler)
	{
		printf("Inverse: Matrix ist singulär\n");
		return 0;
	}
	// Die angehängte Einheitsmatrix Matrix hat sich jetzt in die inverse Matrix umgewandelt
	// Umkopieren auf die Zielmatrix
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv->matrix_entry[i][j] = a->matrix_entry[i][n + j];
		}
	}
	return 1;
}

double myNewAbsAreSoHard(double weakAbs) 
{
	if (weakAbs < 0)
		weakAbs *= -1;
	return weakAbs;
}

int matrix_valid(Matrix* x) 
{
	for (int i = 0; i < x->row_size; i++)
	{
		for (int j = 0; j < x->col_size; j++)
		{
			if (x->matrix_entry[i][j] != x->matrix_entry[i][j])
				return 0;
		}
	}

	return 1;
}

//---------------------------------------------------
//	calculate transpose of matrix
Matrix* matrix_transpose(Matrix* c) 
{
	int i, j;
	Matrix* b = matrix_alloc(c->col_size, c->row_size);

	for (i = 0; i < c->col_size; i++)
		for (j = 0; j < c->row_size; j++)
			b->matrix_entry[i][j] = c->matrix_entry[j][i];

	return b;
}


Matrix* RK4(float sampleTime, Matrix* x_cont, Matrix* x_temp)
{
	Matrix* gammat = matrix_alloc(4, 1);
	Matrix* k1 = matrix_alloc(4, 1);
	Matrix* k2 = matrix_alloc(4, 1);
	Matrix* k3 = matrix_alloc(4, 1);
	Matrix* k4 = matrix_alloc(4, 1);
	Matrix* K = matrix_alloc(4, 4);

	gammat->matrix_entry[0][0] = 1.0f / 6;
	gammat->matrix_entry[1][0] = 1.0f / 3;
	gammat->matrix_entry[2][0] = 1.0f / 3;
	gammat->matrix_entry[3][0] = 1.0f / 6;

	k1->matrix_entry[0][0] = sampleTime * x_cont->matrix_entry[0][0];
	k1->matrix_entry[1][0] = sampleTime * x_cont->matrix_entry[1][0];
	k1->matrix_entry[2][0] = sampleTime * x_cont->matrix_entry[2][0];
	k1->matrix_entry[3][0] = sampleTime * x_cont->matrix_entry[3][0];

	k2->matrix_entry[0][0] = sampleTime * (x_cont->matrix_entry[0][0] + 0.5f * k1->matrix_entry[0][0]);
	k2->matrix_entry[1][0] = sampleTime * (x_cont->matrix_entry[1][0] + 0.5f * k1->matrix_entry[1][0]);
	k2->matrix_entry[2][0] = sampleTime * (x_cont->matrix_entry[2][0] + 0.5f * k1->matrix_entry[2][0]);
	k2->matrix_entry[3][0] = sampleTime * (x_cont->matrix_entry[3][0] + 0.5f * k1->matrix_entry[3][0]);

	k3->matrix_entry[0][0] = sampleTime * (x_cont->matrix_entry[0][0] + 0.5f * k2->matrix_entry[0][0]);
	k3->matrix_entry[1][0] = sampleTime * (x_cont->matrix_entry[1][0] + 0.5f * k2->matrix_entry[1][0]);
	k3->matrix_entry[2][0] = sampleTime * (x_cont->matrix_entry[2][0] + 0.5f * k2->matrix_entry[2][0]);
	k3->matrix_entry[3][0] = sampleTime * (x_cont->matrix_entry[3][0] + 0.5f * k2->matrix_entry[3][0]);

	k4->matrix_entry[0][0] = sampleTime * (x_cont->matrix_entry[0][0] + k3->matrix_entry[0][0]);
	k4->matrix_entry[1][0] = sampleTime * (x_cont->matrix_entry[1][0] + k3->matrix_entry[1][0]);
	k4->matrix_entry[2][0] = sampleTime * (x_cont->matrix_entry[2][0] + k3->matrix_entry[2][0]);
	k4->matrix_entry[3][0] = sampleTime * (x_cont->matrix_entry[3][0] + k3->matrix_entry[3][0]);

	K->matrix_entry[0][0] = k1->matrix_entry[0][0];
	K->matrix_entry[1][0] = k1->matrix_entry[1][0];
	K->matrix_entry[2][0] = k1->matrix_entry[2][0];
	K->matrix_entry[3][0] = k1->matrix_entry[3][0];

	K->matrix_entry[0][1] = k2->matrix_entry[0][0];
	K->matrix_entry[1][1] = k2->matrix_entry[1][0];
	K->matrix_entry[2][1] = k2->matrix_entry[2][0];
	K->matrix_entry[3][1] = k2->matrix_entry[3][0];

	K->matrix_entry[0][2] = k3->matrix_entry[0][0];
	K->matrix_entry[1][2] = k3->matrix_entry[1][0];
	K->matrix_entry[2][2] = k3->matrix_entry[2][0];
	K->matrix_entry[3][2] = k3->matrix_entry[3][0];

	K->matrix_entry[0][3] = k4->matrix_entry[0][0];
	K->matrix_entry[1][3] = k4->matrix_entry[1][0];
	K->matrix_entry[2][3] = k4->matrix_entry[2][0];
	K->matrix_entry[3][3] = k4->matrix_entry[3][0];

	Matrix* K_gammat = matrix_multiply(K, gammat);
	Matrix* x_discrete = matrix_add(x_temp, K_gammat);

	matrix_free(k1);
	matrix_free(k2);
	matrix_free(k3);
	matrix_free(k4);
	matrix_free(K);
	matrix_free(gammat);
	matrix_free(K_gammat);

	return x_discrete;
}


