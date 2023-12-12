#include "fastMult.h"
#include <stdlib.h>

void fastMult3x3MatMat(double *A, double *B, double *C, size_t length) {
  size_t i;

  length*=9;
  for (i=0; i<length; i+=9) {
    C[i  ]=A[i  ]*B[i  ]+A[i+3]*B[i+1]+A[i+6]*B[i+2];
    C[i+1]=A[i+1]*B[i  ]+A[i+4]*B[i+1]+A[i+7]*B[i+2];
    C[i+2]=A[i+2]*B[i  ]+A[i+5]*B[i+1]+A[i+8]*B[i+2];
    C[i+3]=A[i  ]*B[i+3]+A[i+3]*B[i+4]+A[i+6]*B[i+5];
    C[i+4]=A[i+1]*B[i+3]+A[i+4]*B[i+4]+A[i+7]*B[i+5];
    C[i+5]=A[i+2]*B[i+3]+A[i+5]*B[i+4]+A[i+8]*B[i+5];
    C[i+6]=A[i  ]*B[i+6]+A[i+3]*B[i+7]+A[i+6]*B[i+8];
    C[i+7]=A[i+1]*B[i+6]+A[i+4]*B[i+7]+A[i+7]*B[i+8];
    C[i+8]=A[i+2]*B[i+6]+A[i+5]*B[i+7]+A[i+8]*B[i+8];
  }
}

void fastMult3x3MatMat_opt(double *A, double *B, double *C, size_t length) {
  size_t i;
  double *pC=C;
  double *pB=B;
  register double b1,b2,b3;
  
  length*=9;
  for (i=0; i<length; i+=9) {
    b1=*(pB++);
    b2=*(pB++);
    b3=*(pB++);
    *(pC++)=A[i  ]*b1+A[i+3]*b2+A[i+6]*b3;
    *(pC++)=A[i+1]*b1+A[i+4]*b2+A[i+7]*b3;
    *(pC++)=A[i+2]*b1+A[i+5]*b2+A[i+8]*b3;
    b1=*(pB++);
    b2=*(pB++);
    b3=*(pB++);
    *(pC++)=A[i  ]*b1+A[i+3]*b2+A[i+6]*b3;
    *(pC++)=A[i+1]*b1+A[i+4]*b2+A[i+7]*b3;
    *(pC++)=A[i+2]*b1+A[i+5]*b2+A[i+8]*b3;
    b1=*(pB++);
    b2=*(pB++);
    b3=*(pB++);
    *(pC++)=A[i  ]*b1+A[i+3]*b2+A[i+6]*b3;
    *(pC++)=A[i+1]*b1+A[i+4]*b2+A[i+7]*b3;
    *(pC++)=A[i+2]*b1+A[i+5]*b2+A[i+8]*b3;
  }
}

void fastMult3x3HatMat(double *v, double *A, double *B, size_t length) {
  size_t i,j;
  for (i=0, j=0; i<length*9; i+=9, j+=3) {
    B[i  ]= v[j+2]*A[i+1]-v[j+1]*A[i+2];
    B[i+1]=-v[j+2]*A[i  ]+v[j  ]*A[i+2];
    B[i+2]= v[j+1]*A[i  ]-v[j  ]*A[i+1];
    B[i+3]= v[j+2]*A[i+4]-v[j+1]*A[i+5];
    B[i+4]=-v[j+2]*A[i+3]+v[j  ]*A[i+5];
    B[i+5]= v[j+1]*A[i+3]-v[j  ]*A[i+4];
    B[i+6]= v[j+2]*A[i+7]-v[j+1]*A[i+8];
    B[i+7]=-v[j+2]*A[i+6]+v[j  ]*A[i+8];
    B[i+8]= v[j+1]*A[i+6]-v[j  ]*A[i+7];
  }
}

void fastMult3x3MatHat(double *A, double *v, double *B, size_t length) {
  size_t i,j;
  for (i=0, j=0; i<length*9; i+=9, j+=3) {
    B[i  ]=-A[i+3]*v[j+2]+A[i+6]*v[j+1];
    B[i+1]=-A[i+4]*v[j+2]+A[i+7]*v[j+1];
    B[i+2]=-A[i+5]*v[j+2]+A[i+8]*v[j+1];
    B[i+3]=A[i  ]*v[j+2]-A[i+6]*v[j  ];
    B[i+4]=A[i+1]*v[j+2]-A[i+7]*v[j  ];
    B[i+5]=A[i+2]*v[j+2]-A[i+8]*v[j  ];
    B[i+6]=-A[i  ]*v[j+1]+A[i+3]*v[j  ];
    B[i+7]=-A[i+1]*v[j+1]+A[i+4]*v[j  ];
    B[i+8]=-A[i+2]*v[j+1]+A[i+5]*v[j  ];
  }
}

void fastMult3x3MatMatMat(double *A, double *B, double *C, double *D, size_t length) {
  double *E;
  size_t i;

  length*=9;
  E=(double *) malloc(sizeof(double)*length);
  
  for (i=0; i<length; i+=9) {
    E[i  ]=A[i  ]*B[i  ]+A[i+3]*B[i+1]+A[i+6]*B[i+2];
    E[i+1]=A[i+1]*B[i  ]+A[i+4]*B[i+1]+A[i+7]*B[i+2];
    E[i+2]=A[i+2]*B[i  ]+A[i+5]*B[i+1]+A[i+8]*B[i+2];
    E[i+3]=A[i  ]*B[i+3]+A[i+3]*B[i+4]+A[i+6]*B[i+5];
    E[i+4]=A[i+1]*B[i+3]+A[i+4]*B[i+4]+A[i+7]*B[i+5];
    E[i+5]=A[i+2]*B[i+3]+A[i+5]*B[i+4]+A[i+8]*B[i+5];
    E[i+6]=A[i  ]*B[i+6]+A[i+3]*B[i+7]+A[i+6]*B[i+8];
    E[i+7]=A[i+1]*B[i+6]+A[i+4]*B[i+7]+A[i+7]*B[i+8];
    E[i+8]=A[i+2]*B[i+6]+A[i+5]*B[i+7]+A[i+8]*B[i+8];
  }

  for (i=0; i<length; i+=9) {
    D[i  ]=E[i  ]*C[i  ]+E[i+3]*C[i+1]+E[i+6]*C[i+2];
    D[i+1]=E[i+1]*C[i  ]+E[i+4]*C[i+1]+E[i+7]*C[i+2];
    D[i+2]=E[i+2]*C[i  ]+E[i+5]*C[i+1]+E[i+8]*C[i+2];
    D[i+3]=E[i  ]*C[i+3]+E[i+3]*C[i+4]+E[i+6]*C[i+5];
    D[i+4]=E[i+1]*C[i+3]+E[i+4]*C[i+4]+E[i+7]*C[i+5];
    D[i+5]=E[i+2]*C[i+3]+E[i+5]*C[i+4]+E[i+8]*C[i+5];
    D[i+6]=E[i  ]*C[i+6]+E[i+3]*C[i+7]+E[i+6]*C[i+8];
    D[i+7]=E[i+1]*C[i+6]+E[i+4]*C[i+7]+E[i+7]*C[i+8];
    D[i+8]=E[i+2]*C[i+6]+E[i+5]*C[i+7]+E[i+8]*C[i+8];
  }

  free(E);
}

void fastMult3x3MatMatMat_opt(double *A, double *B, double *C, double *D, size_t length) {
  double e1,e2,e3;
  size_t i;

  length*=9;
  for (i=0; i<length; i+=9) {
    e1=A[i  ]*B[i  ]+A[i+3]*B[i+1]+A[i+6]*B[i+2];
    e2=A[i  ]*B[i+3]+A[i+3]*B[i+4]+A[i+6]*B[i+5];
    e3=A[i  ]*B[i+6]+A[i+3]*B[i+7]+A[i+6]*B[i+8];
    D[i  ]=e1*C[i  ]+e2*C[i+1]+e3*C[i+2];
    D[i+3]=e1*C[i+3]+e2*C[i+4]+e3*C[i+5];
    D[i+6]=e1*C[i+6]+e2*C[i+7]+e3*C[i+8];

    e1=A[i+1]*B[i  ]+A[i+4]*B[i+1]+A[i+7]*B[i+2];
    e2=A[i+1]*B[i+3]+A[i+4]*B[i+4]+A[i+7]*B[i+5];
    e3=A[i+1]*B[i+6]+A[i+4]*B[i+7]+A[i+7]*B[i+8];
    D[i+7]=e1*C[i+6]+e2*C[i+7]+e3*C[i+8];
    D[i+4]=e1*C[i+3]+e2*C[i+4]+e3*C[i+5];
    D[i+1]=e1*C[i  ]+e2*C[i+1]+e3*C[i+2];

    e1=A[i+2]*B[i  ]+A[i+5]*B[i+1]+A[i+8]*B[i+2];
    e2=A[i+2]*B[i+3]+A[i+5]*B[i+4]+A[i+8]*B[i+5];
    e3=A[i+2]*B[i+6]+A[i+5]*B[i+7]+A[i+8]*B[i+8];
    D[i+5]=e1*C[i+3]+e2*C[i+4]+e3*C[i+5];
    D[i+2]=e1*C[i  ]+e2*C[i+1]+e3*C[i+2];
    D[i+8]=e1*C[i+6]+e2*C[i+7]+e3*C[i+8];
  }
}

void fastMult3x3HatNormSqMat(double *v, double *A, double *B, size_t length) {
  double c;
  size_t i,j;

  length*=9;
  for (i=0, j=0; i<length; i+=9, j+=3) {
    c=A[i  ]*v[j  ]+A[i+1]*v[j+1]+A[i+2]*v[j+2];
    B[i  ]=c*v[j  ]-A[i  ];
    B[i+1]=c*v[j+1]-A[i+1];
    B[i+2]=c*v[j+2]-A[i+2];
    c=A[i+3]*v[j  ]+A[i+4]*v[j+1]+A[i+5]*v[j+2];
    B[i+3]=c*v[j  ]-A[i+3];
    B[i+4]=c*v[j+1]-A[i+4];
    B[i+5]=c*v[j+2]-A[i+5];
    c=A[i+6]*v[j  ]+A[i+7]*v[j+1]+A[i+8]*v[j+2];
    B[i+6]=c*v[j  ]-A[i+6];
    B[i+7]=c*v[j+1]-A[i+7];
    B[i+8]=c*v[j+2]-A[i+8];
  }
}

void fastMult3x3MatHatNormSq(double *A, double *v, double *B, size_t length) {
  double c;
  size_t i,j;

  length*=9;
  for (i=0, j=0; i<length; i+=9, j+=3) {
    c=A[i  ]*v[j  ]+A[i+3]*v[j+1]+A[i+6]*v[j+2];
    B[i  ]=c*v[j  ]-A[i ];
    B[i+3]=c*v[j+1]-A[i+3];
    B[i+6]=c*v[j+2]-A[i+6];
   
    c=A[i+1]*v[j  ]+A[i+4]*v[j+1]+A[i+7]*v[j+2];
    B[i+1]=c*v[j  ]-A[i+1];
    B[i+4]=c*v[j+1]-A[i+4];
    B[i+7]=c*v[j+2]-A[i+7];

    c=A[i+2]*v[j  ]+A[i+5]*v[j+1]+A[i+8]*v[j+2];
    B[i+2]=c*v[j  ]-A[i+2];
    B[i+5]=c*v[j+1]-A[i+5];
    B[i+8]=c*v[j+2]-A[i+8];
  }
}

void fastMult3x3MatIdxTranspMatMatIdx(double *A, double *eA, double *B, double *C, double *eC, double *D, size_t length) {
  size_t i,k,eeA,eeC;
  double *E;

  length*=9;
  E=(double *) malloc(sizeof(double)*length);

  for (i=0, k=0; i<length; i+=9, ++k) {
    eeA=((size_t) eA[k]-1)*9;
    E[i  ]=A[eeA  ]*B[i  ]+A[eeA+1]*B[i+1]+A[eeA+2]*B[i+2];
    E[i+1]=A[eeA+3]*B[i  ]+A[eeA+4]*B[i+1]+A[eeA+5]*B[i+2];
    E[i+2]=A[eeA+6]*B[i  ]+A[eeA+7]*B[i+1]+A[eeA+8]*B[i+2];
    E[i+3]=A[eeA  ]*B[i+3]+A[eeA+1]*B[i+4]+A[eeA+2]*B[i+5];
    E[i+4]=A[eeA+3]*B[i+3]+A[eeA+4]*B[i+4]+A[eeA+5]*B[i+5];
    E[i+5]=A[eeA+6]*B[i+3]+A[eeA+7]*B[i+4]+A[eeA+8]*B[i+5];
    E[i+6]=A[eeA  ]*B[i+6]+A[eeA+1]*B[i+7]+A[eeA+2]*B[i+8];
    E[i+7]=A[eeA+3]*B[i+6]+A[eeA+4]*B[i+7]+A[eeA+5]*B[i+8];
    E[i+8]=A[eeA+6]*B[i+6]+A[eeA+7]*B[i+7]+A[eeA+8]*B[i+8];
  }
  
  for (i=0, k=0; i<length; i+=9, ++k) {
    eeC=((size_t) eC[k]-1)*9;
    D[i  ]=E[i  ]*C[eeC  ]+E[i+3]*C[eeC+1]+E[i+6]*C[eeC+2];
    D[i+1]=E[i+1]*C[eeC  ]+E[i+4]*C[eeC+1]+E[i+7]*C[eeC+2];
    D[i+2]=E[i+2]*C[eeC  ]+E[i+5]*C[eeC+1]+E[i+8]*C[eeC+2];
    D[i+3]=E[i  ]*C[eeC+3]+E[i+3]*C[eeC+4]+E[i+6]*C[eeC+5];
    D[i+4]=E[i+1]*C[eeC+3]+E[i+4]*C[eeC+4]+E[i+7]*C[eeC+5];
    D[i+5]=E[i+2]*C[eeC+3]+E[i+5]*C[eeC+4]+E[i+8]*C[eeC+5];
    D[i+6]=E[i  ]*C[eeC+6]+E[i+3]*C[eeC+7]+E[i+6]*C[eeC+8];
    D[i+7]=E[i+1]*C[eeC+6]+E[i+4]*C[eeC+7]+E[i+7]*C[eeC+8];
    D[i+8]=E[i+2]*C[eeC+6]+E[i+5]*C[eeC+7]+E[i+8]*C[eeC+8];
  }

  free(E);
}

void fastMult3x3MatIdxTranspMatMatIdx_opt(double *A, double *eA, double *B, double *C, double *eC, double *D, size_t length) {
  size_t i,k,eeA,eeC;
  double E1,E2,E3;

  length*=9;
  for (i=0, k=0; i<length; i+=9, ++k) {
    eeA=((size_t) eA[k]-1)*9;
    eeC=((size_t) eC[k]-1)*9;
    E1=A[eeA  ]*B[i  ]+A[eeA+1]*B[i+1]+A[eeA+2]*B[i+2];
    E2=A[eeA  ]*B[i+3]+A[eeA+1]*B[i+4]+A[eeA+2]*B[i+5];
    E3=A[eeA  ]*B[i+6]+A[eeA+1]*B[i+7]+A[eeA+2]*B[i+8];
    D[i  ]=E1*C[eeC  ]+E2*C[eeC+1]+E3*C[eeC+2];
    D[i+3]=E1*C[eeC+3]+E2*C[eeC+4]+E3*C[eeC+5];
    D[i+6]=E1*C[eeC+6]+E2*C[eeC+7]+E3*C[eeC+8];

    E1=A[eeA+3]*B[i  ]+A[eeA+4]*B[i+1]+A[eeA+5]*B[i+2];
    E2=A[eeA+3]*B[i+3]+A[eeA+4]*B[i+4]+A[eeA+5]*B[i+5];
    E3=A[eeA+3]*B[i+6]+A[eeA+4]*B[i+7]+A[eeA+5]*B[i+8];
    D[i+1]=E1*C[eeC  ]+E2*C[eeC+1]+E3*C[eeC+2];
    D[i+4]=E1*C[eeC+3]+E2*C[eeC+4]+E3*C[eeC+5];
    D[i+7]=E1*C[eeC+6]+E2*C[eeC+7]+E3*C[eeC+8];

    E1=A[eeA+6]*B[i  ]+A[eeA+7]*B[i+1]+A[eeA+8]*B[i+2];
    E2=A[eeA+6]*B[i+3]+A[eeA+7]*B[i+4]+A[eeA+8]*B[i+5];
    E3=A[eeA+6]*B[i+6]+A[eeA+7]*B[i+7]+A[eeA+8]*B[i+8];
    D[i+2]=E1*C[eeC  ]+E2*C[eeC+1]+E3*C[eeC+2];
    D[i+5]=E1*C[eeC+3]+E2*C[eeC+4]+E3*C[eeC+5];
    D[i+8]=E1*C[eeC+6]+E2*C[eeC+7]+E3*C[eeC+8];
  }
}
