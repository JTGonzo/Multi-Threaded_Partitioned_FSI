/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "StVenantKirchhoffMaterial.h"

/*************************************************************************/
void StVenantKirchhoffMaterial_forces(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    
    double* myRrows    = mxGetPr(plhs[0]);
    double* myRcoef    = mxGetPr(plhs[1]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double* elements  = mxGetPr(prhs[4]);
        
    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }
    
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myRrows,myRcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        
        double gradphi[NumQuadPoints][dim][nln];

        double GradV[dim][dim];
        double GradUh[NumQuadPoints][dim][dim];
        
        double F[NumQuadPoints][dim][dim];
        double E[NumQuadPoints][dim][dim];
        double P_Uh[dim][dim];
        double FE[NumQuadPoints][dim][dim];
                
        double traceE[NumQuadPoints];
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
            
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[q][d1][d2] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[q][d1][d2] = GradUh[q][d1][d2] + U_h[e_k] * gradphi[q][d2][k];
                    }
                    F[q][d1][d2]  = Id[d1][d2] + GradUh[q][d1][d2];
                }
            }

            MatrixProductAlphaT1(dim, 1.0, F[q], F[q], E[q] );
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    E[q][d1][d2] = 0.5 * ( E[q][d1][d2] - Id[d1][d2] );
                }
            }
            traceE[q] = Trace(dim, E[q]);
            
            MatrixProduct(dim, F[q], E[q], FE[q] );

        }

        int ii = 0;
        int a, i_c, j_c;
        
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }
                
                double rloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[i_c][d2] = gradphi[q][d2][a];
                    }
                    
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            P_Uh[d1][d2] =  ( 2 * mu * FE[q][d1][d2] + lambda * traceE[q] * F[q][d1][d2] );
                        }
                    }
                    
                    rloc  = rloc + Mdot( dim, GradV, P_Uh) * w[q];
                }
                                            
                myRrows[ie*nln*dim+ii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                myRcoef[ie*nln*dim+ii] = rloc*detjac[ie];
                ii = ii + 1;
            }
        }
    }
}

/*************************************************************************/
void StVenantKirchhoffMaterial_jacobianFast2D(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double* elements  = mxGetPr(prhs[4]);
        
    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        
        double GradV[dim][dim];
        double GradU[dim][dim];
        double GradUh[NumQuadPoints][dim][dim];
        double F[NumQuadPoints][dim][dim];
        double E[NumQuadPoints][dim][dim];
        double traceE[NumQuadPoints];

        double gradphi[NumQuadPoints][dim][nln];
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
            
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[q][d1][d2] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[q][d1][d2] = GradUh[q][d1][d2] + U_h[e_k] * gradphi[q][d2][k];
                    }
                    F[q][d1][d2]  = Id[d1][d2] + GradUh[q][d1][d2];
                }
            }
            
            MatrixProductAlphaT1(dim, 1.0, F[q], F[q], E[q] );
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    E[q][d1][d2] = 0.5 * ( E[q][d1][d2] - Id[d1][d2] );
                }
            }
            traceE[q] = Trace(dim, E[q]);
        }

        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;
        
        double aloc[nln][dim][nln][dim];
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 2; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 2; j_c = j_c + 1 )
                    {
                        aloc[a][i_c][b][j_c]  = 0.0;
                    }
                }
            }
        }
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* loop over test functions --> a */
            for (a = 0; a < nln; a = a + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    
                    aloc[a][0][b][0] += ( gradphi[q][0][a]*(F[q][0][0]*(lambda*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b]) + 2.0*F[q][0][0]*mu*gradphi[q][0][b]) + gradphi[q][0][b]*(2.0*E[q][0][0]*mu + lambda*traceE[q]) + 2.0*F[q][0][1]*mu*((F[q][0][0]*gradphi[q][1][b])/2.0 + (F[q][0][1]*gradphi[q][0][b])/2.0) + 2.0*E[q][1][0]*mu*gradphi[q][1][b]) + gradphi[q][1][a]*(F[q][0][1]*(lambda*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b]) + 2.0*F[q][0][1]*mu*gradphi[q][1][b]) + gradphi[q][1][b]*(2.0*E[q][1][1]*mu + lambda*traceE[q]) + 2.0*F[q][0][0]*mu*((F[q][0][0]*gradphi[q][1][b])/2.0 + (F[q][0][1]*gradphi[q][0][b])/2.0) + 2.0*E[q][0][1]*mu*gradphi[q][0][b]) ) * w[q];
                    
                    aloc[a][0][b][1] += ( gradphi[q][0][a]*(F[q][0][0]*(lambda*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b]) + 2.0*F[q][1][0]*mu*gradphi[q][0][b]) + 2.0*F[q][0][1]*mu*((F[q][1][0]*gradphi[q][1][b])/2.0 + (F[q][1][1]*gradphi[q][0][b])/2.0)) + gradphi[q][1][a]*(F[q][0][1]*(lambda*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b]) + 2.0*F[q][1][1]*mu*gradphi[q][1][b]) + 2.0*F[q][0][0]*mu*((F[q][1][0]*gradphi[q][1][b])/2.0 + (F[q][1][1]*gradphi[q][0][b])/2.0)) ) * w[q];
                    
                    aloc[a][1][b][0] += ( gradphi[q][0][a]*(F[q][1][0]*(lambda*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b]) + 2.0*F[q][0][0]*mu*gradphi[q][0][b]) + 2.0*F[q][1][1]*mu*((F[q][0][0]*gradphi[q][1][b])/2.0 + (F[q][0][1]*gradphi[q][0][b])/2.0)) + gradphi[q][1][a]*(F[q][1][1]*(lambda*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b]) + 2.0*F[q][0][1]*mu*gradphi[q][1][b]) + 2.0*F[q][1][0]*mu*((F[q][0][0]*gradphi[q][1][b])/2.0 + (F[q][0][1]*gradphi[q][0][b])/2.0)) ) * w[q];
                    
                    aloc[a][1][b][1] += ( gradphi[q][0][a]*(F[q][1][0]*(lambda*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b]) + 2.0*F[q][1][0]*mu*gradphi[q][0][b]) + gradphi[q][0][b]*(2.0*E[q][0][0]*mu + lambda*traceE[q]) + 2.0*F[q][1][1]*mu*((F[q][1][0]*gradphi[q][1][b])/2.0 + (F[q][1][1]*gradphi[q][0][b])/2.0) + 2.0*E[q][1][0]*mu*gradphi[q][1][b]) + gradphi[q][1][a]*(F[q][1][1]*(lambda*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b]) + 2.0*F[q][1][1]*mu*gradphi[q][1][b]) + gradphi[q][1][b]*(2.0*E[q][1][1]*mu + lambda*traceE[q]) + 2.0*F[q][1][0]*mu*((F[q][1][0]*gradphi[q][1][b])/2.0 + (F[q][1][1]*gradphi[q][0][b])/2.0) + 2.0*E[q][0][1]*mu*gradphi[q][0][b]) ) * w[q];

                }
            }
        }

        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 2; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 2; j_c = j_c + 1 )
                    {
                        myArows[ie*nln2*4+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*4+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*4+iii] = aloc[a][i_c][b][j_c]*detjac[ie];
                        iii = iii + 1;
                    }
                }
            }
        }
    }
        
}