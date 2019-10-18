#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>
#include <complex>

#include "Real.h"
#include "Complex.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "LinearContext.h"
#include "Utility.h"

namespace DDG
{
   extern thread_local LinearContext context;

   const int maxEigIter = 20;
   // number of iterations used to solve eigenvalue problems

   template <class T>
   SparseMatrix<T> :: SparseMatrix( size_t m_, size_t n_ )
   // initialize an mxn matrix
   : m( m_ ),
     n( n_ ),
     cData( NULL )
   {}

   template <class T>
   SparseMatrix<T> :: SparseMatrix( const SparseMatrix<T>& B )
   // copy constructor
   : cData( NULL )
   {
      *this = B;
   }

   template <class T>
   SparseMatrix<T> :: ~SparseMatrix( void )
   // destructor
   {
      if( cData )
      {
         cholmod_l_free_sparse( &cData, context );
      }
   }

   template <class T>
   const SparseMatrix<T>& SparseMatrix<T> :: operator=( const SparseMatrix<T>& B )
   // copies B
   {
      if( cData )
      {
         cholmod_l_free_sparse( &cData, context );
         cData = NULL;
      }

      m = B.m;
      n = B.n;
      data = B.data;

      return *this;
   }

   template <class T>
   SparseMatrix<T> SparseMatrix<T> :: transpose( void ) const
   {
      SparseMatrix<T> AT( n, m );

      for( const_iterator e  = begin();
                          e != end();
                          e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         T Aij = e->second;

         AT(j,i) = Aij.bar();
      }

      return AT;
   }

   template <class T>
   SparseMatrix<T> SparseMatrix<T> :: operator*( const SparseMatrix<T>& B ) const
   // returns product of this matrix with sparse B
   {
      const SparseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nColumns() == B.nRows() );

      // collect nonzeros in each row
      std::vector< std::vector< int > > Bcol( B.nRows() );
      std::vector< std::vector<  T  > > Bval( B.nRows() );
      for( const_iterator e  = B.begin();
                          e != B.end();
                          e ++ )
      {
         int row = e->first.second;
         int col = e->first.first;
         T val = e->second;

         Bcol[ row ].push_back( col );
         Bval[ row ].push_back( val );
      }

      // multiply C = A*B
      SparseMatrix<T> C( A.nRows(), B.nColumns() );
      for( const_iterator e  = begin();
                          e != end();
                          e ++ )
      {
         int i = e->first.second;
         int j = e->first.first;

         for( size_t nn = 0; nn < Bcol[j].size(); nn++ )
         {
            int k = Bcol[j][nn];

            C( i, k ) += e->second * Bval[j][nn];
         }
      }

      return C;
   }

   template <class T>
   DenseMatrix<T> SparseMatrix<T> :: operator*( const DenseMatrix<T>& B ) const
   // returns product of this matrix with dense B
   {
      const SparseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nColumns() == B.nRows() );

      // multiply C = A*B
      DenseMatrix<T> C( A.nRows(), B.nColumns() );
      for( const_iterator e  = begin();
                          e != end();
                          e ++ )
      {
         int i = e->first.second;
         int j = e->first.first;

         for( int k = 0; k < B.nColumns(); k++ )
         {
            C( i, k ) += e->second * B( j, k );
         }
      }

      return C;
   }

   template <class T>
   void SparseMatrix<T> :: operator*=( const T& c )
   {
      for( iterator e  = begin();
                    e != end();
                    e++ )
      {
         e->second *= c;
      }
   }

   template <class T>
   void SparseMatrix<T> :: operator/=( const T& c )
   {
      for( iterator e  = begin();
                    e != end();
                    e++ )
      {
         e->second /= c;
      }
   }

   template <class T>
   void SparseMatrix<T> :: operator+=( const SparseMatrix<T>& B )
   // adds B to this matrix
   {
      SparseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nRows() == B.nRows() );
      assert( A.nColumns() == B.nColumns() );

      for( const_iterator e  = B.begin();
                          e != B.end();
                          e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         const T& Bij( e->second );

         A( i, j ) += Bij;
      }
   }

   template <class T>
   void SparseMatrix<T> :: operator-=( const SparseMatrix<T>& B )
   // subtracts B from this matrix
   {
      SparseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nRows() == B.nRows() );
      assert( A.nColumns() == B.nColumns() );

      for( const_iterator e  = B.begin();
                          e != B.end();
                          e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         const T& Bij( e->second );

         A( i, j ) -= Bij;
      }
   }

   template <class T>
   SparseMatrix<T> SparseMatrix<T> :: operator+( const SparseMatrix<T>& B ) const
   // returns sum of this matrix with B
   {
      SparseMatrix<T> C( nRows(), nColumns() );

      C += *this;
      C += B;

      return C;
   }

   template <class T>
   SparseMatrix<T> SparseMatrix<T> :: operator-( const SparseMatrix<T>& B ) const
   // returns sum of this matrix with B
   {
      SparseMatrix<T> C( nRows(), nColumns() );

      C += *this;
      C -= B;

      return C;
   }

   template <class T>
   SparseMatrix<T> operator*( const T& c, const SparseMatrix<T>& A )
   {
      SparseMatrix<T> cA = A;

      for( typename SparseMatrix<T>::iterator e  = cA.begin();
                                              e != cA.end();
                                              e++ )
      {
         e->second = c * e->second;
      }

      return cA;
   }

   template <class T>
   SparseMatrix<T> operator*( const SparseMatrix<T>& A, const T& c )
   {
      SparseMatrix<T> Ac = A;

      Ac *= c;

      return Ac;
   }

   template <class T>
   SparseMatrix<T> operator/( const SparseMatrix<T>& A, T c )
   {
      SparseMatrix<T> Ac = A;

      Ac /= c;

      return Ac;
   }

   template <class T>
   void SparseMatrix<T>::resize( size_t m_, size_t n_)
   {
      m = m_;
      n = n_;

      data.clear();
   }

   template <class T>
   size_t SparseMatrix<T> :: nRows( void ) const
   // returns the number of rows
   {
      return m;
   }

   template <class T>
   size_t SparseMatrix<T> :: nColumns( void ) const
   // returns the number of columns
   {
      return n;
   }

   template <class T>
   int SparseMatrix<T> :: length( void ) const
   // returns the size of the largest dimension
   {
      return std::max( m, n );
   }

   template <class T>
   void SparseMatrix<T> :: zero( const T& val )
   // sets all nonzero elements val
   {
      for( iterator i  = begin();
                    i != end();
                    i ++ )
      {
         i->second = val;
      }
   }

   template <class T>
   void SparseMatrix<T> :: invertDiagonal( void )
   {
      SparseMatrix<T>& A( *this );

      for( int i = 0; i < std::max( m, n ); i++ )
      {
         A( i, i ) = A( i, i ).inv();
      }
   }

   template <class T>
   SparseMatrix<T> SparseMatrix<T> :: identity( int N )
   {
      SparseMatrix<T> I( N, N );

      for( int i = 0; i < N; i++ )
      {
         I( i, i ) = 1.;
      }

      return I;
   }

   template <class T>
   DenseMatrix<T> SparseMatrix<T> :: full( void ) const
   // converts to a dense matrix
   {
      const int maxSize = 1048576;
      if( m*n > maxSize )
      {
         std::cerr << "Error: refusing to convert sparse to dense (too big!)" << "\n";
         exit( 1 );
      }

      const SparseMatrix<T>& A( *this );
      DenseMatrix<T> B( m, n );

      for( int i = 0; i < m; i++ )
      for( int j = 0; j < n; j++ )
      {
         B( i, j ) = A( i, j );
      }

      return B;
   }

   template <class T>
   cholmod_sparse* SparseMatrix<T> :: to_cholmod( void )
   {
      if( cData )
      {
         cholmod_l_free_sparse( &cData, context );
         cData = NULL;
      }

      allocateSparse();

      // build compressed matrix (note that EntryMap stores entries in column-major order)
       double* pr =  (double*) cData->x;
      SuiteSparse_long* ir = (SuiteSparse_long*) cData->i;
      SuiteSparse_long* jc = (SuiteSparse_long*) cData->p;
      int i = 0;
      int j = -1;
      for( const_iterator e  = begin();
                          e != end();
                          e ++ )
      {
         int c = e->first.first;
         if( c != j )
         {
            for( int k = j+1; k <= c; k++ )
            {
               jc[k] = i;
            }
            j = c;
         }

         ir[i] = e->first.second;
         setEntry( e, i, pr );

         i++;
      }
      for( int k = j+1; k <= n; k++ )
      {
         jc[k] = i;
      }

      return cData;
   }

   template <>
   cholmod_sparse* SparseMatrix<Quaternion> :: to_cholmod( void );

   template <class T>
   T& SparseMatrix<T> :: operator()( int row, int col )
   {
      EntryIndex index( col, row );
      const_iterator entry = data.find( index );

      if( entry == end())
      {
         data[ index ] = T( 0. );
      }

      return data[ index ];
   }

   template <class T>
   T SparseMatrix<T> :: operator()( int row, int col ) const
   {
      EntryIndex index( col, row );
      const_iterator entry = data.find( index );

      if( entry == end())
      {
         return T( 0. );
      }

      return entry->second;
   }

   template <class T>
   typename SparseMatrix<T>::iterator SparseMatrix<T> :: begin( void )
   {
      return data.begin();
   }

   template <class T>
   typename SparseMatrix<T>::const_iterator SparseMatrix<T> :: begin( void ) const
   {
      return data.begin();
   }

   template <class T>
   typename SparseMatrix<T>::iterator SparseMatrix<T> :: end( void )
   {
      return data.end();
   }

   template <class T>
   typename SparseMatrix<T>::const_iterator SparseMatrix<T> :: end( void ) const
   {
      return data.end();
   }

   template <class T>
   void SparseMatrix<T> :: shift( double c )
   // adds c times the identity matrix to this matrix
   {
      assert( m == n );
      SparseMatrix<T>& A( *this );

      for( int i = 0; i < m; i++ )
      {
         A( i, i ) += c;
      }
   }

   template <class T>
   void solvePositiveDefinite( SparseMatrix<T>& A,
                                DenseMatrix<T>& x,
                                DenseMatrix<T>& b )
   // solves the positive definite sparse linear system Ax = b using sparse Cholesky factorization
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
      cholmod_sparse* Ac = A.to_cholmod();
      Ac->stype = 1;
      cholmod_factor* L = cholmod_l_analyze( Ac, context );
      cholmod_l_factorize( Ac, L, context );
      x = cholmod_l_solve( CHOLMOD_A, L, b.to_cholmod(), context );

      if( L ) cholmod_l_free_factor( &L, context );
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[chol] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[chol] max residual: " << residual( A, x, b ) << "\n";
#endif
   }

   template <class T>
   void backsolvePositiveDefinite(  SparseFactor<T>& L,
                                     DenseMatrix<T>& x,
                                     DenseMatrix<T>& b )
   // backsolves the prefactored positive definite sparse linear system LL'x = b
   {
      x = cholmod_l_solve( CHOLMOD_A, L.to_cholmod(), b.to_cholmod(), context );
   }

   template <class T>
   void smallestEig( SparseMatrix<T>& A,
                      DenseMatrix<T>& x,
                      bool ignoreConstantVector )
   // solves A x = lambda x for the smallest nonzero eigenvalue lambda
   // A must be symmetric; x is used as an initial guess
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
      for( int iter = 0; iter < maxEigIter; iter++ )
      {
         solve( A, x, x );
         if( ignoreConstantVector )
         {
            x.removeMean();
         }
         x.normalize();
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[eig] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[eig] max residual: " << residual( A, x ) << "\n";
#endif
   }

   template <class T>
   void smallestEig( SparseMatrix<T>& A,
                     SparseMatrix<T>& B,
                      DenseMatrix<T>& x )
   // solves A x = lambda B x for the smallest nonzero generalized eigenvalue lambda
   // A and B must be symmetric; x is used as an initial guess
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif

      DenseMatrix<T> y;
      std::mt19937 gen;
      x.randomize(gen);

      for( int iter = 0; iter < maxEigIter; iter++ )
      {
         y = B*x;
         solve( A, x, y );
         x /= sqrt( inner( B*x, x ).re );
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[eig] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[eig] max residual: " << residual( A, B, x ) << "\n";
#endif
   }

   template <class T>
   void smallestEigPositiveDefinite( SparseMatrix<T>& A,
                                      DenseMatrix<T>& x,
                                      bool ignoreConstantVector )
   // solves A x = lambda x for the smallest nonzero eigenvalue lambda
   // A must be positive (semi-)definite; x is used as an initial guess
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
      SparseFactor<T> L;
      L.build( A );

      for( int iter = 0; iter < maxEigIter; iter++ )
      {
         backsolvePositiveDefinite( L, x, x );
         if( ignoreConstantVector )
         {
            x.removeMean();
         }
         x.normalize();
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[eig] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[eig] max residual: " << residual( A, x ) << "\n";
#endif
   }

   template <class T>
   void smallestEigPositiveDefinite( SparseMatrix<T>& A,
                                     SparseMatrix<T>& B,
                                      DenseMatrix<T>& x )
   // solves A x = lambda x for the smallest nonzero eigenvalue lambda
   // A must be positive (semi-)definite; x is used as an initial guess
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif

      SparseFactor<T> L;
      L.build( A );

      std::mt19937 gen;
      x.randomize(gen);

      for( int iter = 0; iter < maxEigIter; iter++ )
      {
         x = B*x;
         backsolvePositiveDefinite( L, x, x );
         x /= sqrt( inner( B*x, x ).norm() );
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[eig] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[eig] max residual: " << residual( A, B, x ) << "\n";
      cout << "[eig] energy: " << inner( A*x, x ) << "\n";
#endif
   }

   template <class T>
   void smallestEigPositiveDefinite( SparseMatrix<T>& A,
                                     SparseMatrix<T>& B,
                                      DenseMatrix<T>& E,
                                      DenseMatrix<T>& x )
   // solves A x = lambda B x for the smallest nonzero eigenvalue lambda orthogonal
   // to the vectors spanned by the columns of E, which must be orthonormal. A must
   // be positive (semi-)definite, B must be symmetric; x is used as an initial guess
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif

      SparseFactor<T> L;
      L.build( A );

      std::mt19937 gen;
      x.randomize(gen);

      for( int iter = 0; iter < maxEigIter; iter++ )
      {
         x = B*x;
         backsolvePositiveDefinite( L, x, x );

         // project out components of E
         for( int k = 0; k < E.nColumns(); k++ )
         {
            double xDotEk = 0.;
            for( int i = 0; i < x.nRows(); i++ )
            {
               xDotEk += x(i) * E(i,k);
            }
            for( int i = 0; i < x.nRows(); i++ )
            {
               x(i) -= xDotEk * E(i,k);
            }
         }

         x /= sqrt( inner( B*x, x ).norm() );
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[eig] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[eig] max residual: " << residual( A, B, x ) << "\n";
      cout << "[eig] energy: " << inner( A*x, x ) << "\n";
#endif
   }

   template <class T>
   void nSmallestEigsPositiveDefinite( SparseMatrix<T>& A,
                                       SparseMatrix<T>& B,
                                       std::vector<DenseMatrix<T> >& V,
                                       std::vector<double>& D,
                                       int nEigs )
   // computes the n smallest eigenvectors and eigenvalues of A with
   // respect to the inner product B; V stores the eigenvectors as
   // columns and D is a list of eigenvalues (in the same order); n
   // is the requested number of eigenpairs
   {
      int nRows = A.nRows();
      V.resize( nEigs );
      D.resize( nEigs );

      SparseFactor<T> L;
      L.build( A );

      std::mt19937 gen;
      for( int k = 0; k < nEigs; k++ )
      {
         V[k] = DenseMatrix<T>( nRows );
         V[k].randomize(gen);

         for( int iter = 0; iter < maxEigIter; iter++ )
         {
            V[k] = B*V[k];
            backsolvePositiveDefinite( L, V[k], V[k] );

            // project out previous vectors
            for( int j = 0; j < k; j++ )
            {
               Real Cjk = inner( V[j], V[k] );
               V[k] -= Cjk * V[j];
            }

            V[k] /= sqrt( inner( B*V[k], V[k] ).norm() );
         }

         D[k] = inner( A*V[k], V[k] );
      }
   }

   template <class T>
   double residual( const SparseMatrix<T>& A,
                    const  DenseMatrix<T>& x,
                    const  DenseMatrix<T>& b )
   // returns the max residual of the linear problem A x = b relative to the largest entry of the solution
   {
      return ( A*x - b ).norm() / b.norm();
   }

   template <class T>
   double residual( const SparseMatrix<T>& A,
                    const  DenseMatrix<T>& x )
   // returns the max residual of the eigenvalue problem A x = lambda x relative to the largest entry of the solution
   {
      T lambda = rayleighQuotient( A, x );
      return (A*x-lambda*x).norm() / x.norm();
   }

   template <class T>
   double residual( const SparseMatrix<T>& A,
                    const SparseMatrix<T>& B,
                    const  DenseMatrix<T>& x )
   // returns the max residual of the generalized eigenvalue problem A x = lambda x relative to the largest entry of the solution
   {
      T lambda = rayleighQuotient( A, B, x );
      return (A*x-lambda*(B*x)).norm() / x.norm();
   }

   template <class T>
   double residual( const SparseMatrix<T>& A,
                    const SparseMatrix<T>& B,
                    const  DenseMatrix<T>& E,
                    const  DenseMatrix<T>& x )
   // returns the max residual of the generalized eigenvalue problem A x = lambda (B - EE^T) x relative to the largest entry of the solution
   {
      T lambda = rayleighQuotient( A, B, E, x );
      return (A*x-lambda*(B*x-E*(E.transpose()*x))).norm() / x.norm();
   }

   template <class T>
   T rayleighQuotient( const SparseMatrix<T>& A,
                       const  DenseMatrix<T>& x )
   // returns <x,Ax>/<x,x>
   {
      return (x.transpose()*(A*x))(0) * (x.transpose()*x)(0).inv();
   }

   template <class T>
   T rayleighQuotient( const SparseMatrix<T>& A,
                       const SparseMatrix<T>& B,
                       const  DenseMatrix<T>& x )
   // returns <Ax,x>/<Bx,x>
   {
      return (x.transpose()*(A*x))(0) * (x.transpose()*(B*x))(0).inv();
   }

   template <class T>
   T rayleighQuotient( const SparseMatrix<T>& A,
                       const SparseMatrix<T>& B,
                       const  DenseMatrix<T>& E,
                       const  DenseMatrix<T>& x )
   // returns <Ax,x>/<(B-EE^T)x,x>
   {
      return (x.transpose()*(A*x))(0) * (x.transpose()*(B*x-E*(E.transpose()*x)))(0).inv();
   }

   template <class T>
   std::ostream& operator<<( std::ostream& os, const SparseMatrix<T>& o)
   {
      os.precision( 3 );

      for( typename SparseMatrix<T>::const_iterator e  = o.begin();
                                                    e != o.end();
                                                    e ++ )
      {
         int row = e->first.second;
         int col = e->first.first;

         os << "( " << row << ", " << col << " ): " << e->second << "\n";
      }

      return os;
   }

   template <class T>
   SparseFactor<T> :: SparseFactor( void )
   : L( NULL )
   {}

   template <class T>
   SparseFactor<T> :: ~SparseFactor( void )
   {
      if( L )
      {
         cholmod_l_free_factor( &L, context );
      }
   }

   template <class T>
   void SparseFactor<T> :: build( SparseMatrix<T>& A )
   {
      if( L )
      {
         cholmod_l_free_factor( &L, context );
         L = NULL;
      }

      cholmod_sparse* Ac = A.to_cholmod();
      Ac->stype = 1;

#ifdef SP_DEBUG
      int t0, t1;
      t0 = clock();
#endif
      L = cholmod_l_analyze( Ac, context );
#ifdef SP_DEBUG
      t1 = clock();
      std::cerr << "analyze: " << seconds(t0,t1) << "s" << std::endl;
#endif

#ifdef SP_DEBUG
      t0 = clock();
#endif
      cholmod_l_factorize( Ac, L, context );
#ifdef SP_DEBUG
      t1 = clock();
      std::cerr << "factorize: " << seconds(t0,t1) << "s" << std::endl;
#endif
   }

   template <class T>
   bool SparseFactor<T> :: valid( void ) const
   {
      if( L == NULL )
      {
         return false;
      }
      return true;
   }

   template <class T>
   cholmod_factor* SparseFactor<T> :: to_cholmod( void )
   {
      return L;
   }
}

