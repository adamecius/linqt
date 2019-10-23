

#ifndef ALGEBRA_FUNCTION
#define ALGEBRA_FUNCTION


//VECTOR FUNCTIONS
void CreateRandomVector( const int dim, const int i0, complex* x );

// VECTOR-VECTOR FUNCTIONS
void axpy(const int dim , complex a,const complex* x, complex* y );

void copy(const int dim , const complex* x, complex* y );

complex dot(const int dim , const complex* x, complex* y ) ;


struct CSRMatrix
{

	void Multiply(const complex a,const complex* x,const complex  b, complex* y );

	integer Dim(){ return dim_;};

	bool ReadCSRMatrix( const std::string input);

	void Optimize( const int numMul )
	{
		mkl_sparse_set_mv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_set_dotmv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_optimize (A);
	}




	void Rescale( const complex a  )
	{
		for(int i = 0; i < values.size(); i++) 
			values[i]=values[i]*a;
	}

	private:
		struct matrix_descr descr;
		sparse_matrix_t A;
//		char descr[6];		
		integer dim_;
		std::vector<complex> values;
		std::vector<integer> columns, rowIndex;
};







void axpy(const int dim , complex a,const complex* x, complex* y )
{
	cblas_zaxpy(dim, &a,x,1,y,1);
}

void copy(const int dim , const complex* x, complex* y )
{
	cblas_zcopy(dim,x,1,y,1);
}

complex dot(const int dim , const complex* x, complex* y )
{
	complex dotc;
	cblas_zdotc_sub (dim ,x,1,y,1, &dotc);
	return dotc;
}


void CreateRandomVector( const int dim, const int i0,const int di, complex* x )
{
	for(int  i = i0 ; i-i0 < dim; i+=di)
	if( i < dim )
		x[i] = exp( complex(0.0,2.0*M_PI*(double)rand()/(double)RAND_MAX ) )/sqrt(dim);
};


void CSRMatrix::Multiply(const complex a,const complex* x,const complex  b, complex* y )
{

	if( mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,
						a, A ,descr,
						x,b,y) != SPARSE_STATUS_SUCCESS )
	std::cout<<"PROBLEMS"<<std::endl;
//	else
//	std::cout<<"SUCESS"<<std::endl;

//	mkl_zcsrmv(&descr[4] , &dim_ ,&dim_ ,&a , &descr[0],&values[0] , &columns[0] , &rowIndex[0] , &rowIndex[1] , x , &b , y );


};

bool CSRMatrix::ReadCSRMatrix( const std::string input)
{
	//OPEN MATRIX FILE
	std::ifstream matrix_file ( input.c_str()) ;

	//READ DIMENSION OF THE MATRIX
	int nnz, dim;
	matrix_file>>dim>>nnz; dim_=dim;

	//CREATE ARRAYS TO STORE THE MATRIX
	values  = std::vector<complex>(nnz);
	columns = std::vector<integer>(nnz);
	rowIndex= std::vector<integer>(dim+1);

	//READ VALUES
	double rev,imv;
	for(int i=0;i < nnz; i++)
	{
		matrix_file>>rev>>imv;
		values[i]= complex(rev,imv);
	}

	//READ COLUMNS
	int col;
	for(int i=0;i < nnz; i++)
	{
		matrix_file>>col;
		columns[i]= col;
	}

		//READ ROW_INDEX_ARRAY
	int rowIdx;
	for(int i=0;i < dim+1; i++)
	{
		matrix_file>>rowIdx;
		rowIndex[i]= rowIdx;
	}
	matrix_file.close();

	//TRANSFORM THE ARRAYS IN A MKL_SPARSE_MATRIX
	descr.type = SPARSE_MATRIX_TYPE_HERMITIAN ;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	if( mkl_sparse_z_create_csr(&A, SPARSE_INDEX_BASE_ZERO,Dim(),Dim(), &rowIndex[0],&rowIndex[1],&columns[0],&values[0]) != SPARSE_STATUS_SUCCESS )
		std::cout<<"PROBLEMS"<<std::endl;
//		else
//		std::cout<<"SUCESS"<<std::endl;
//	descr[0] = 'H'; descr[1] = 'U'; descr[2] = 'N'; 
//	descr[3] = 'C'; descr[4] = 'N'; 

	

	};



#endif
