#include "SparseMatrix.h"

namespace DDG
{
   template <>
   const SparseMatrix<Real>& SparseMatrix<Real> :: operator=( cholmod_sparse* B )
   {
      assert( B );
      assert( B->xtype == CHOLMOD_REAL );

      if( cData )
      {
         cholmod_l_free_sparse( &cData, context );
      }
      cData = B;

      m = cData->nrow;
      n = cData->ncol;
      resize( m, n );

      double* pr = (double*) cData->x;
      SuiteSparse_long* ir = (SuiteSparse_long*) cData->i;
      SuiteSparse_long* jc = (SuiteSparse_long*) cData->p;

      // iterate over columns
      for( int col = 0; col < n; col++ )
      {
         // iterate over nonzero rows
         for( auto k = jc[col]; k < jc[col+1]; k++ )
         {
            int row = (int) ir[k];

            (*this)( row, col ) = pr[k];
         }
      }

      return *this;
   }

   template <>
   const SparseMatrix<Complex>& SparseMatrix<Complex> :: operator=( cholmod_sparse* B )
   {
      assert( B );
      assert( B->xtype == CHOLMOD_COMPLEX );

      if( cData )
      {
         cholmod_l_free_sparse( &cData, context );
      }
      cData = B;

      m = cData->nrow;
      n = cData->ncol;
      resize( m, n );

      double* pr = (double*) cData->x;
      SuiteSparse_long* ir = (SuiteSparse_long*) cData->i;
      SuiteSparse_long* jc = (SuiteSparse_long*) cData->p;

      // iterate over columns
      for( int col = 0; col < n; col++ )
      {
         // iterate over nonzero rows
         for( auto k = jc[col]; k < jc[col+1]; k++ )
         {
            int row = (int) ir[k];

            (*this)( row, col ) = Complex( pr[k*2+0], pr[k*2+1] );
         }
      }

      return *this;
   }

   template <>
   cholmod_sparse* SparseMatrix<Quaternion> :: to_cholmod( void )
   {
      SparseMatrix<Real> A( m*4, n*4 );

      for( const_iterator e  = begin();
                          e != end();
                          e ++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         const Quaternion& q( e->second );

         A(i*4+0,j*4+0) =  q[0]; A(i*4+0,j*4+1) = -q[1]; A(i*4+0,j*4+2) = -q[2]; A(i*4+0,j*4+3) = -q[3];
         A(i*4+1,j*4+0) =  q[1]; A(i*4+1,j*4+1) =  q[0]; A(i*4+1,j*4+2) = -q[3]; A(i*4+1,j*4+3) =  q[2];
         A(i*4+2,j*4+0) =  q[2]; A(i*4+2,j*4+1) =  q[3]; A(i*4+2,j*4+2) =  q[0]; A(i*4+2,j*4+3) = -q[1];
         A(i*4+3,j*4+0) =  q[3]; A(i*4+3,j*4+1) = -q[2]; A(i*4+3,j*4+2) =  q[1]; A(i*4+3,j*4+3) =  q[0];
      }

      if( cData != NULL )
      {
         cholmod_l_free_sparse( &cData, context );
      }
      cData = cholmod_l_copy_sparse( A.to_cholmod(), context );
      return cData;
   }

   template <>
   void SparseMatrix<Real> :: allocateSparse( void )
   {
      auto nzmax = data.size();
      int sorted = true;
      int packed = true;
      int stype = 0;
      cData = cholmod_l_allocate_sparse( m, n, nzmax, sorted, packed, stype, CHOLMOD_REAL, context );
   }

   template <>
   void SparseMatrix<Complex> :: allocateSparse( void )
   {
      auto nzmax = data.size();
      int sorted = true;
      int packed = true;
      int stype = 0;
      cData = cholmod_l_allocate_sparse( m, n, nzmax, sorted, packed, stype, CHOLMOD_COMPLEX, context );
   }

   template <>
   void SparseMatrix<Quaternion> :: allocateSparse( void )
   {
      auto nzmax = data.size();
      int sorted = true;
      int packed = true;
      int stype = 0;
      cData = cholmod_l_allocate_sparse( m*4, n*4, nzmax, sorted, packed, stype, CHOLMOD_REAL, context );
   }

   template <>
   void SparseMatrix<Real> :: setEntry( const_iterator e, int i, double* pr )
   {
      pr[i] = e->second;
   }

   template <>
   void SparseMatrix<Complex> :: setEntry( const_iterator e, int i, double* pr )
   {
      pr[i*2+0] = e->second.re;
      pr[i*2+1] = e->second.im;
   }
}
