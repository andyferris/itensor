#ifndef __lapack_wrap_h
#define __lapack_wrap_h

//
// Headers and typedefs
//
#ifdef PLATFORM_macos

#include <Accelerate/Accelerate.h>
typedef __CLPK_integer
LAPACK_INT;
typedef __CLPK_doublereal
LAPACK_REAL;
typedef __CLPK_doublecomplex
LAPACK_COMPLEX;

#elif PLATFORM_acml

#include "acml.h"
typedef int
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef doublecomplex
LAPACK_COMPLEX;

#elif PLATFORM_mkl

#include "mkl_lapack.h"
typedef MKL_INT
LAPACK_INT;
typedef double
LAPACK_REAL;
typedef MKL_Complex16
LAPACK_COMPLEX;

#endif


//
// dsyev
//
void inline
dsyev_wrapper(char* jobz,        //if jobz=='V', compute eigs and evecs
              char* uplo,        //if uplo=='U', read from upper triangle of A
              LAPACK_INT* n,     //number of cols of A
              LAPACK_REAL* A,    //symmetric matrix A
              LAPACK_INT* lda,   //size of A (usually same as n)
              LAPACK_REAL* eigs, //eigenvalues on return
              LAPACK_INT* info)  //error info
    {
    LAPACK_INT lwork = max(1,3*(*n)-1);
    LAPACK_REAL work[lwork];
#ifdef PLATFORM_acml
    dsyev_(jobz,uplo,n,A,lda,eigs,work,&lwork,info,1,1);
#else
    dsyev_(jobz,uplo,n,A,lda,eigs,work,&lwork,info);
#endif
    }

void inline
zgesdd_wrapper(char *jobz,           //char* specifying how much of U, V to compute
                                     //choosing *jobz=='S' computes min(m,n) cols of U, V
               LAPACK_INT *m,        //number of rows of input matrix *A
               LAPACK_INT *n,        //number of cols of input matrix *A
               LAPACK_COMPLEX *A,    //contents of input matrix A
               LAPACK_REAL *s,       //on return, singular values of A
               LAPACK_COMPLEX *u,    //on return, unitary matrix U
               LAPACK_COMPLEX *vt,   //on return, unitary matrix V transpose
               LAPACK_INT *info)
    {
    LAPACK_INT l = min(*m,*n),
               g = max(*m,*n);
    LAPACK_INT lwork = l*l+2*l+g+100;
    LAPACK_COMPLEX work[lwork];
    LAPACK_REAL rwork[5*l*(1+l)];
    LAPACK_INT iwork[8*l];
    zgesdd_(jobz,m,n,A,m,s,u,m,vt,n,work,&lwork,rwork,iwork,info);
    }

//
// dgeqrf
//
// QR factorization of a real matrix A
//
void inline
dgeqrf_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_REAL* A,    //matrix A
                                  //on return upper triangle contains R
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_REAL* tau,  //scalar factors of elementary reflectors
                                  //length should be min(m,n)
               LAPACK_INT* info)  //error info
    {
    int lwork = max(1,4*max(*n,*m));
    LAPACK_REAL work[lwork]; 
    dgeqrf_(m,n,A,lda,tau,work,&lwork,info);
    }

//
// dorgqr
//
// Generates Q from output of QR factorization routine dgeqrf (see above)
//
void inline
dorgqr_wrapper(LAPACK_INT* m,     //number of rows of A
               LAPACK_INT* n,     //number of cols of A
               LAPACK_INT* k,     //number of elementary reflectors, typically min(m,n)
               LAPACK_REAL* A,    //matrix A, as returned from "A" argument of dgeqrf
                                  //on return contains Q
               LAPACK_INT* lda,   //size of A (usually same as n)
               LAPACK_REAL* tau,  //scalar factors as returned by dgeqrf
               LAPACK_INT* info)  //error info
    {
    int lwork = max(1,4*max(*n,*m));
    LAPACK_REAL work[lwork]; 
    dorgqr_(m,n,k,A,lda,tau,work,&lwork,info);
    }

//
// zheev
//
// Eigenvalues and eigenvectors of complex Hermitian matrix A
//
void inline
zheev_wrapper(char* jobz,           //if 'V', compute both eigs and evecs
                                    //if 'N', only eigenvalues
              char* uplo,           //if 'U', use upper triangle of A
              LAPACK_INT* n,        //number of cols of A
              LAPACK_COMPLEX* A,    //matrix A, on return contains eigenvectors
              LAPACK_INT* lda,      //size of A (usually same as n)
              LAPACK_REAL* d,       //eigenvalues on return
              LAPACK_COMPLEX* work, //complex workspace array
              LAPACK_INT* lwork,    //size of work
              LAPACK_REAL* rwork,   //real workspace array
              LAPACK_INT* info)  //error info
    {
    zheev_(jobz,uplo,n,A,lda,d,work,lwork,rwork,info);
    }

//
// dsygv
//
// Eigenvalues and eigenvectors of generalized eigenvalue problem
// A*x = lambda*B*x
// A and B must be symmetric
// B must be positive definite
//
void inline
dsygv_wrapper(char* jobz,           //if 'V', compute both eigs and evecs
                                    //if 'N', only eigenvalues
              char* uplo,           //if 'U', use upper triangle of A
              LAPACK_INT* n,        //number of cols of A
              LAPACK_REAL* A,       //matrix A, on return contains eigenvectors
              LAPACK_REAL* B,       //matrix B
              LAPACK_REAL* d,       //eigenvalues on return
              LAPACK_INT* info)  //error info
    {
    int itype = 1;
    LAPACK_INT lwork = max(1,3*(*n)-1);//max(1, 1+6*N+2*N*N);
    LAPACK_REAL work[lwork];
    dsygv_(&itype,jobz,uplo,n,A,n,B,n,d,work,&lwork,info);
    }

//
// dgeev
//
// Eigenvalues and eigenvectors of real, square matrix A
// A can be a general real matrix, not assumed symmetric
//
void inline
dgeev_wrapper(char* jobvl,          //if 'V', compute left eigenvectors, else 'N'
              char* jobvr,          //if 'V', compute right eigenvectors, else 'N'
              LAPACK_INT* n,        //number of rows/cols of A
              LAPACK_REAL* A,       //matrix A, on return contains eigenvectors
              LAPACK_REAL* dr,      //real parts of eigenvalues
              LAPACK_REAL* di,      //imaginary parts of eigenvalues
              LAPACK_REAL* vl,      //left eigenvectors on return
              LAPACK_REAL* vr,      //right eigenvectors on return
              LAPACK_INT* info)  //error info
    {
    int nevecl = (*jobvl == 'V' ? *n : 1);
    int nevecr = (*jobvr == 'V' ? *n : 1);
    LAPACK_INT lwork = max(1,4*(*n));
    LAPACK_REAL work[lwork];
    dgeev_(jobvl,jobvr,n,A,n,dr,di,vl,&nevecl,vr,&nevecr,work,&lwork,info);
    }

#endif
