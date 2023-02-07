#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>


/****************************************************/
/* convienence method for evaluating an R function  */
/* that returns 1,2 or 3 and takes in two rows of a */
/* matrix.                                          */
/* @param Rfun  the function                        */
/* @param Rvect1 the first vector for the comparison*/
/* @param Rvect2 the second vector for th comparison*/
/**/
/* @return the return value of Rfun                 */
/****************************************************/
double run_fun(SEXP Rfun, SEXP Rvect1, SEXP Rvect2) {
  SEXP e, result;
  
  PROTECT(e = lang3(Rfun, Rvect1, Rvect2));
  result = eval(e, R_GlobalEnv);
  UNPROTECT(1);
  
  return (REAL(result)[0]);
}


/*******************************************************************/
/* convienence method to extract a row from a matrix as an R vector*/
/*                                                                 */
/* @param matrix the matrix to extract from                        */
/* @param row    the row to extract                                */
/* @return an SEXP pointing to an R vector                         */
/*******************************************************************/
SEXP extract_row(SEXP matrix, int row) {
  int i;
  SEXP rc;
  SEXP dim = getAttrib(matrix, R_DimSymbol);
  int rows = INTEGER(dim)[0];
  int cols = INTEGER(dim)[1];
  
  PROTECT(rc = allocVector(REALSXP, cols));
  for (i=0; i<cols;i++) {
    REAL(rc)[i] = REAL(matrix)[row+i*rows];
  }
  UNPROTECT(1);
  return rc;
}

/***********************************************************************/
/* experimental pi function. Can take a generic function in and        */
/* calculates the pi function based on the return of fun               */
/*                                                                     */
/* @param Rpostmat the matrix with the data in it.                     */
/* @param Rfun the function to evaluate the relation between points    */
/* @param Rr the maximum distances to look at                          */
/* @param Rr_min the minimum distances                                 */
/* @param Rinds  indices into the original array, to help with bootstrapping*/
/* @param Rxcol the column containing the x coordinate                 */
/* @param Rycol the column containing the y coordinate                 */
/* @param Renvr the enviroment for function evaluatoin                 */
/***********************************************************************/
SEXP get_pi (SEXP Rpostmat,
             SEXP Rfun,
             SEXP Rr,
             SEXP Rr_low,
             SEXP Rinds,
             SEXP Rxcol,
             SEXP Rycol) {
  
  SEXP rc = R_NilValue;
  int i,j,k;
  double dist;
  int num_cnt, denom_cnt; /*counters for those filling conditions*/
int f_ans; /*used to hold the result of the function*/


/*turn all of the R stuff passed in to the type of stuff we can
referene in C*/
double *r = REAL(Rr);
double *r_low = REAL(Rr_low);
int *inds = INTEGER(Rinds);
int xcol = INTEGER(Rxcol)[0]-1;
int ycol = INTEGER(Rycol)[0]-1;

SEXP postmat_dim = getAttrib(Rpostmat, R_DimSymbol);
double *postmat = REAL(Rpostmat);
int rows = INTEGER(postmat_dim)[0];

/*some sanity checking*/
if (!isFunction(Rfun)) error("Rfun must be a function");

/*prepare the return information*/
PROTECT(rc=allocVector(REALSXP, length(Rr)));

/*repeat calculation for all r*/
for (i=0;i<length(Rr);i++) {
  //Rprintf("%1.1f,%1.1f\n", r_low[i],r[i]); //DEBUG
  /*zero out counts*/
  num_cnt = 0;
  denom_cnt = 0;
  
  /*might be faster to have some scraning criteria, but
  we will loop throufh every pair...j is from */
  for (j=0;j<rows;j++) {
    
    /*k is to*/
    for (k=0; k<rows;k++) {
      /*do not compare someone with themself*/
      if (inds[k] == inds[j]) continue;
      
      /*calculate the distance*/
      dist = sqrt(pow(postmat[j+xcol*rows] - postmat[k+xcol*rows],2) +
      pow(postmat[j+ycol*rows] - postmat[k+ycol*rows],2));
      
      if ((dist>r[i]) | (dist<r_low[i])) continue;
      
      /*call the user supplied function*/
      f_ans = (int)run_fun(Rfun,
               extract_row(Rpostmat,j),
               extract_row(Rpostmat,k));
      
      /*update the counts appropriately*/
      if (f_ans==1) {
        denom_cnt++;
        num_cnt++;
      } else if (f_ans==2) {
        denom_cnt++;
      }
    }
  }
  //Rprintf("%d/%d\n",num_cnt,denom_cnt); // DEBUG
  REAL(rc)[i] = (double)num_cnt/denom_cnt;
}

UNPROTECT(1);

return(rc);

}

/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_typed (int *type,
                   double *x,
                   double *y,
                   int *len,
                   int *typeA,
                   int *typeB,
                   double *r_low,
                   double *r,
                   int *len_r,
                   int *inds,
                   double *rc) {
  
  int i,j,k;
  int num_cnt, denom_cnt; /*counters for those filling conditions*/
double dist;

/*repeat the calculation for all r*/
for (i=0;i<*len_r;i++) {
  //Rprintf("%1.1f,%1.1f\n", r_low[i],r[i]); //DEBUG
  /*zero out counts*/
  num_cnt = 0;
  denom_cnt = 0;
  
  if (*typeA != -1) {
    
    for (j=0;j<*len;j++) {
      
      if (type[j] != *typeA) continue;
      
      for (k=0;k<*len;k++) {
        /*ignore pairs of the same type*/
        if (inds[k]==inds[j]) continue;
        
        dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
        if ((dist<=r[i])  & (dist>=r_low[i])) denom_cnt++;
        
        if (type[k] != *typeB) continue;
        if ((dist<=r[i])  & (dist>=r_low[i])) num_cnt++;
      }
    }
    
  } else {
    Rprintf("To be implemented\n");
    return;
  }
  
  //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
  rc[i] = (double)num_cnt/denom_cnt;
}
}



/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_typed_fast (int *type,
                        double *x,
                        double *y,
                        int *len,
                        int *typeA,
                        int *typeB,
                        double *r_low,
                        double *r,
                        int *len_r,
                        int *inds,
                        double *rc) {
  
  int i,j,k;
  int num_cnt, denom_cnt; /*counters for those filling conditions*/
double dist;

/*repeat the calculation for all r*/
for (i=0;i<*len_r;i++) {
  //Rprintf("%1.1f,%1.1f\n", r_low[i],r[i]); //DEBUG
  /*zero out counts*/
  num_cnt = 0;
  denom_cnt = 0;
  
  if (*typeA != -1) {
    
    for (j=0;j<*len;j++) {
      
      if (type[j] != *typeA) continue;
      
      for (k=0;k<*len;k++) {
        /*ignore pairs of the same type*/
        if (inds[k]==inds[j]) continue;
        
        dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
        if ((dist<=r[i])  & (dist>=r_low[i])) denom_cnt++;
        
        if (type[k] != *typeB) continue;
        if ((dist<=r[i])  & (dist>=r_low[i])) num_cnt++;
      }
    }
    
  } else {
    Rprintf("To be implemented\n");
    return;
  }
  
  //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
  rc[i] = (double)num_cnt/denom_cnt;
}
}




/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_typed_wts (int *type,
                       double *x,
                       double *y,
                       double *weight,
                       int *len,
                       int *typeA,
                       int *typeB,
                       double *r_low,
                       double *r,
                       int *len_r,
                       int *inds,
                       double *rc) {
  
  int i,j,k;
  int num_cnt, denom_cnt; /*counters for those filling conditions*/
double dist;

/*repeat the calculation for all r*/
for (i=0;i<*len_r;i++) {
  //Rprintf("%1.1f,%1.1f\n", r_low[i],r[i]); //DEBUG
  /*zero out counts*/
  num_cnt = 0;
  denom_cnt = 0;
  
  if (*typeA != -1) {
    
    for (j=0;j<*len;j++) {
      if (type[j] != *typeA) continue;
      
      for (k=0;k<*len;k++) {
        /*ignore pairs of the same type*/
        if (inds[k]==inds[j]) continue;
        
        dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
        if ((dist<=r[i])  & (dist>=r_low[i])) denom_cnt = denom_cnt + weight[k]*weight[j];
        
        if (type[k] != *typeB) continue;
        if ((dist<=r[i])  & (dist>=r_low[i])) num_cnt = num_cnt + weight[k]*weight[j];
      }
    }
    
  } else {
    Rprintf("To be implemented\n");
    return;
  }
  //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
  rc[i] = (double)num_cnt/denom_cnt;
}
}



/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_clustsurvey (double *p,
                         double *x,
                         double *y,
                         int *s,
                         int *len,
                         double *r_low,
                         double *r,
                         int *len_r,
                         int *inds,
                         double *rc) {
  
int i,j,k;
int skip_r, skip_phi;  /*counters for those filling conditions*/
double dist;
double si, sj, sum_sj;
double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, phi, sum_phi;
double mean_probsuscept;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
  
  sum_phi = 0;
  skip_r = 1;
  sum_typeAi = 0;
  si = 0;
  
  for (i=0;i<*len;i++) {
    
    skip_phi = 1;
    phi = 0;
    sum_sj = 0;		
    num_typeAi = 0;
    sum_typeAj = 0;			
    
    // Skip cluster if no individuals with the characteristic
    if (p[i] == 0) continue;
    
    // Cycle through paired clusters with cluster j
    //  - need to add s and r to each cluster within the correct distance to the sums
    //  - need to sum the sum_typeAj; 
    for (j=0;j<*len;j++) {
      num_typeAj = 0;
      sj = 0;
      
      dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      
      if ((dist<=r[k])  & (dist>=r_low[k])) {
        sj = s[j];  // number of individuals in the cluster j
        num_typeAj = p[j]*sj;  // number of individuals of type A in cluster i (type A = susceptible)
        
        // Subtract 1 from total number in the cluster and 1 from number susceptible if looking at contact within cluster.
        // Cannot infect self;
        if (inds[i]==inds[j]) {
          sj = s[j] - 1;
          num_typeAj = num_typeAj - 1;
        }
        
        sum_typeAj = sum_typeAj + num_typeAj;
        sum_sj = sum_sj + sj;  // sum of s[j] at distance r from cluster i
        skip_phi = 0;
      }
    }
    
    if (skip_phi == 1) continue;
    
    si = s[i];  // number of individuals in the cluster i
    num_typeAi = p[i]*si;  // number of individuals of type A in cluster i (type A = susceptibles)
    sum_typeAi = sum_typeAi + num_typeAi;	  // sum of number of type A in primary clusters [i]
    mean_probsuscept = sum_typeAj / sum_sj;
    
    phi = num_typeAi*mean_probsuscept;
    sum_phi = sum_phi + phi;
    skip_r = 0;
  }
  
  //Skip if no phi values for that distance of r
  if (skip_r == 1) {
    rc[k] = NA_REAL;
  } else {
    rc[k] = (double)sum_phi/sum_typeAi;
  }
  
  
  rc[k] = (double)sum_phi/sum_typeAi;
}
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_clustsurvey_wts (double *p,
                             double *x,
                             double *y,
                             int *s,
                             double *weight,
                             int *len,
                             double *r_low,
                             double *r,
                             int *len_r,
                             int *inds,
                             double *rc) {
  
int i,j,k;
int skip_r, skip_phi;  /*counters for those filling conditions*/
double dist;
double si, sj, sum_sj;
double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, phi, sum_phi;
double mean_probsuscept;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
  
  sum_phi = 0;
  skip_r = 1;
  sum_typeAi = 0;
  si = 0;
  
  for (i=0;i<*len;i++) {
    
    skip_phi = 1;
    phi = 0;
    sum_sj = 0;		
    num_typeAi = 0;
    sum_typeAj = 0;			
    
    // Skip cluster if no individuals with the characteristic
    if (p[i] == 0) continue;
    
    // Cycle through paired clusters with cluster j
    //  - need to add s and r to each cluster within the correct distance to the sums
    //  - need to sum the sum_typeAj; 
    for (j=0;j<*len;j++) {
      num_typeAj = 0;
      sj = 0;
      
      dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      
      if ((dist<=r[k])  & (dist>=r_low[k])) {
        sj = s[j]*weight[j];  // weighted number of individuals in the cluster j
        num_typeAj = p[j]*sj;  // weighted number of individuals of type A in cluster i (type A = susceptible)
        
        // Subtract 1 from total number in the cluster and 1 from number susceptible if looking at contact within cluster.
        // Cannot infect self;
        if (inds[i]==inds[j]) {
          sj = s[j] - 1;
          num_typeAj = num_typeAj - 1;
        }
        
        sum_typeAj = sum_typeAj + num_typeAj;
        sum_sj = sum_sj + sj;  // sum of weighted s[j] at distance r from cluster i
        skip_phi = 0;
      }
    }
    
    if (skip_phi == 1) continue;
    
    si = s[i]*weight[i];  // weighted number of individuals in the cluster i
    num_typeAi = p[i]*si;  // weighted number of individuals of type A in cluster i (type A = susceptibles)
    sum_typeAi = sum_typeAi + num_typeAi;	  // sum of weighted number of type A in primary clusters [i]
    mean_probsuscept = sum_typeAj / sum_sj;
    
    phi = num_typeAi*mean_probsuscept;
    sum_phi = sum_phi + phi;
    skip_r = 0;
  }
  
  //Skip if no phi values for that distance of r
  if (skip_r == 1) {
    rc[k] = NA_REAL;
  } else {
    rc[k] = (double)sum_phi/sum_typeAi;
  }
  
  rc[k] = (double)sum_phi/sum_typeAi;
}
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_clustsurvey_wts2 (double *p,
                              double *x,
                              double *y,
                              int *s,
                              double *weight,
                              int *len,
                              double *r_low,
                              double *r,
                              int *len_r,
                              int *inds,
                              double *rc) {
  
  int i,j,k;
  int skip_phi;  /*counters for those filling conditions*/
//int skip_r, skip_phi;  /*counters for those filling conditions*/
int count_clust;
double dist;
double si, sum_susc, n_susc, mean_nsusc;
double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, sum_weights, phi, sum_phi;
double mean_probsuscept;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
  
  sum_phi = 0;
  //skip_r = 1;
  sum_typeAi = 0;
  si = 0;
  
  for (i=0;i<*len;i++) {
    
    skip_phi = 1;
    phi = 0;
    sum_susc = 0;		
    num_typeAi = 0;
    sum_typeAj = 0;			
    sum_weights = 0;
    count_clust = 0;
    mean_nsusc = 0;
    
    // Skip cluster if no individuals with the characteristic
    if (p[i] == 0) continue;
    
    // Cycle through paired clusters with cluster j
    //  - need to add s and r to each cluster within the correct distance to the sums
    //  - need to sum the sum_typeAj; 
    for (j=0;j<*len;j++) {
      num_typeAj = 0;
      n_susc = 0;
      
      //ignore pairs of the same cluster
      //if (inds[i]==inds[j]) continue;  //Don't ignore
      
      dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      
      if ((dist<=r[k])  & (dist>=r_low[k])) {
        num_typeAj = p[j]*s[j]*weight[j];  // weighted number of individuals of type A in cluster i (type A = susceptible)
        n_susc = p[j]*s[j];
        sum_typeAj = sum_typeAj + num_typeAj;
        sum_weights = sum_weights + weight[j];
        sum_susc = sum_susc + n_susc;
        skip_phi = 0;
        count_clust = count_clust + 1;
      }
    }
    
    if (skip_phi == 1) continue;
    
    si = s[i]*weight[i];  // weighted number of individuals in the cluster i
    num_typeAi = p[i]*si;  // weighted number of individuals of type A in cluster i (type A = susceptibles)
    sum_typeAi = sum_typeAi + num_typeAi;	  // sum of weighted number of type A in primary clusters [i]
    mean_nsusc = sum_susc / count_clust;
      mean_probsuscept = sum_typeAj / (sum_weights * mean_nsusc);
    
    phi = num_typeAi * mean_probsuscept;
    sum_phi = sum_phi + phi;
    //skip_r = 0;
  }
  
  //Skip if no phi values for that distance of r
  /*		if (skip_r == 1) {
  rc[k] = NA_REAL;
} else {
  rc[k] = (double)sum_phi/sum_typeAi;
}
  */
  
  rc[k] = (double)sum_phi/sum_typeAi;
}
}








/*-----------------------------------------------------------------------------------*/





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_typed_clustsurvey (
                         int *type,
                         double *p,
                         double *x,
                         double *y,
                         int *s,
                         int *len,
                         int *typeA,
                         int *typeB,
                         double *r_low,
                         double *r,
                         int *len_r,
                         int *inds,
                         double *rc) {
    
int i,j,k;
int skip_r, skip_phi;  /*counters for those filling conditions*/
double dist;
double si, sj, sum_sj;
double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, phi, sum_phi;
double mean_probsuscept;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
    
    sum_phi = 0;
    skip_r = 1;
    sum_typeAi = 0;
    si = 0;
    
    // Loop through all clusters as primary clusters
    for (i=0;i<*len;i++) {
        
        skip_phi = 1;
        phi = 0;
        sum_sj = 0;		
        num_typeAi = 0;
        sum_typeAj = 0;			
        
        // Skip cluster if no individuals with the characteristic
        if (p[i] == 0) continue;
        if (type[i] != *typeA) continue;
        
        // Cycle through paired clusters with cluster j
        //  - need to add s and r to each cluster within the correct distance to the sums
        //  - need to sum the sum_typeAj; 
        for (j=0;j<*len;j++) {
            num_typeAj = 0;
            sj = 0;
            
            dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
            
            if ((dist<=r[k])  & (dist>=r_low[k])) {
                sj = s[j];  // number of individuals in the cluster j
                num_typeAj = p[j]*sj;  // number of individuals of type A in cluster i (type A = susceptible)
                
                // Subtract 1 from total number in the cluster and 1 from number susceptible if looking at contact within cluster.
                // Cannot infect self;
                if (inds[i]==inds[j]) {
                    sj = s[j] - 1;
                    num_typeAj = num_typeAj - 1;
                }
                
                sum_typeAj = sum_typeAj + num_typeAj;
                sum_sj = sum_sj + sj;  // sum of s[j] at distance r from cluster i
                skip_phi = 0;
            }
        }
        
        if (skip_phi == 1) continue;
        
        si = s[i];  // number of individuals in the cluster i
        num_typeAi = p[i]*si;  // number of individuals of type A in cluster i (type A = susceptibles)
        sum_typeAi = sum_typeAi + num_typeAi;	  // sum of number of type A in primary clusters [i]
        mean_probsuscept = sum_typeAj / sum_sj;
        
        phi = num_typeAi*mean_probsuscept;
        sum_phi = sum_phi + phi;
        skip_r = 0;
    }
    
    //Skip if no phi values for that distance of r
    if (skip_r == 1) {
        rc[k] = NA_REAL;
    } else {
        rc[k] = (double)sum_phi/sum_typeAi;
    }
    
    
    rc[k] = (double)sum_phi/sum_typeAi;
}
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_typed_clustsurvey_wts (double *p,
                             int *type,
                             double *x,
                             double *y,
                             int *s,
                             double *weight,
                             int *len,
                             int *typeA,
                             int *typeB,
                             double *r_low,
                             double *r,
                             int *len_r,
                             int *inds,
                             double *rc) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
double dist;
double si, sj, sum_sj;
double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, phi, sum_phi;
double mean_probsuscept;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
    
    sum_phi = 0;
    skip_r = 1;
    sum_typeAi = 0;
    si = 0;
    
    for (i=0;i<*len;i++) {
        
        skip_phi = 1;
        phi = 0;
        sum_sj = 0;		
        num_typeAi = 0;
        sum_typeAj = 0;			
        
        // Skip cluster if no individuals with the characteristic
        if (p[i] == 0) continue;
        if (type[i] != *typeA) continue;
        
        // Cycle through paired clusters with cluster j
        //  - need to add s and r to each cluster within the correct distance to the sums
        //  - need to sum the sum_typeAj; 
        for (j=0;j<*len;j++) {
            num_typeAj = 0;
            sj = 0;
            
            dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
            
            if ((dist<=r[k])  & (dist>=r_low[k])) {
                sj = s[j]*weight[j];  // weighted number of individuals in the cluster j
                num_typeAj = p[j]*sj;  // weighted number of individuals of type A in cluster i (type A = susceptible)
                
                // Subtract 1 from total number in the cluster and 1 from number susceptible if looking at contact within cluster.
                // Cannot infect self;
                if (inds[i]==inds[j]) {
                    sj = s[j] - 1;
                    num_typeAj = num_typeAj - 1;
                }
                
                sum_typeAj = sum_typeAj + num_typeAj;
                sum_sj = sum_sj + sj;  // sum of weighted s[j] at distance r from cluster i
                skip_phi = 0;
            }
        }
        
        if (skip_phi == 1) continue;
        
        si = s[i]*weight[i];  // weighted number of individuals in the cluster i
        num_typeAi = p[i]*si;  // weighted number of individuals of type A in cluster i (type A = susceptibles)
        sum_typeAi = sum_typeAi + num_typeAi;	  // sum of weighted number of type A in primary clusters [i]
        mean_probsuscept = sum_typeAj / sum_sj;
        
        phi = num_typeAi*mean_probsuscept;
        sum_phi = sum_phi + phi;
        skip_r = 0;
    }
    
    //Skip if no phi values for that distance of r
    if (skip_r == 1) {
        rc[k] = NA_REAL;
    } else {
        rc[k] = (double)sum_phi/sum_typeAi;
    }
    
    rc[k] = (double)sum_phi/sum_typeAi;
}
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_clustsurvey_wts2 (int *type,
                              double *p,
                              double *x,
                              double *y,
                              int *s,
                              double *weight,
                              int *len,
                              int *typeA,
                              int *typeB,
                              double *r_low,
                              double *r,
                              int *len_r,
                              int *inds,
                              double *rc) {
    
    int i,j,k;
    int skip_phi;  /*counters for those filling conditions*/
//int skip_r, skip_phi;  /*counters for those filling conditions*/
int count_clust;
double dist;
double si, sum_susc, n_susc, mean_nsusc;
double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, sum_weights, phi, sum_phi;
double mean_probsuscept;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
    
    sum_phi = 0;
    //skip_r = 1;
    sum_typeAi = 0;
    si = 0;
    
    for (i=0;i<*len;i++) {
        
        skip_phi = 1;
        phi = 0;
        sum_susc = 0;		
        num_typeAi = 0;
        sum_typeAj = 0;			
        sum_weights = 0;
        count_clust = 0;
        mean_nsusc = 0;
        
        // Skip cluster if no individuals with the characteristic
        if (p[i] == 0) continue;
        if (type[i] != *typeA) continue;
        
        // Cycle through paired clusters with cluster j
        //  - need to add s and r to each cluster within the correct distance to the sums
        //  - need to sum the sum_typeAj; 
        for (j=0;j<*len;j++) {
            num_typeAj = 0;
            n_susc = 0;
            
            //ignore pairs of the same cluster
            //if (inds[i]==inds[j]) continue;  //Don't ignore
            
            dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
            
            if ((dist<=r[k])  & (dist>=r_low[k])) {
                num_typeAj = p[j]*s[j]*weight[j];  // weighted number of individuals of type A in cluster i (type A = susceptible)
                n_susc = p[j]*s[j];
                sum_typeAj = sum_typeAj + num_typeAj;
                sum_weights = sum_weights + weight[j];
                sum_susc = sum_susc + n_susc;
                skip_phi = 0;
                count_clust = count_clust + 1;
            }
        }
        
        if (skip_phi == 1) continue;
        
        si = s[i]*weight[i];  // weighted number of individuals in the cluster i
        num_typeAi = p[i]*si;  // weighted number of individuals of type A in cluster i (type A = susceptibles)
        sum_typeAi = sum_typeAi + num_typeAi;	  // sum of weighted number of type A in primary clusters [i]
        mean_nsusc = sum_susc / count_clust;
            mean_probsuscept = sum_typeAj / (sum_weights * mean_nsusc);
        
        phi = num_typeAi * mean_probsuscept;
        sum_phi = sum_phi + phi;
        //skip_r = 0;
    }
    
    //Skip if no phi values for that distance of r
    /*		if (skip_r == 1) {
    rc[k] = NA_REAL;
} else {
    rc[k] = (double)sum_phi/sum_typeAi;
}
    */
    
    rc[k] = (double)sum_phi/sum_typeAi;
}
}























/*****************************************************************/
/*tau function for generic functions                             */
/**/
/*****************************************************************/
SEXP get_tau (SEXP Rpostmat,
              SEXP Rfun,
              SEXP Rr,
              SEXP Rr_low,
              SEXP Rinds,
              SEXP Rxcol,
              SEXP Rycol) {
  
  int i;
  SEXP divisor;
  SEXP rc;
  SEXP highR;
  SEXP lowR;
  
  PROTECT(highR=allocVector(REALSXP, 1));
  PROTECT(lowR=allocVector(REALSXP, 1));
  REAL(highR)[0] = R_PosInf;
  REAL(lowR)[0] = 0;
  PROTECT(divisor=get_pi(Rpostmat,Rfun, highR, lowR, Rinds, Rxcol, Rycol));
  PROTECT(rc=get_pi(Rpostmat, Rfun, Rr, Rr_low, Rinds, Rxcol, Rycol));
  for (i=0; i<length(Rr);i++) {
    REAL(rc)[i] = REAL(rc)[i]/REAL(divisor)[0];
  }
  UNPROTECT(4);
  
  return(rc);
}


/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_typed (int *type,
                    double *x,
                    double *y,
                    int *len,
                    int *typeA,
                    int *typeB,
                    double *r_low,
                    double *r,
                    int *len_r,
                    int *inds,
                    double *rc) {
  
  int i = 0;
  double divisor;
  double tmp_r_low = 0;
  double tmp_r = DBL_MAX;
  int tmp_len_r = 1;
  
  /*get the divisor in the pi function*/
  get_pi_typed(type,x,y,len,typeA,typeB,&tmp_r_low,&tmp_r,
               &tmp_len_r,inds,&divisor);
  
  /*get the main pi function*/
  get_pi_typed(type,x,y,len,typeA,typeB,r_low,r,len_r,inds,rc);
  
  for (i = 0; i < *len_r; i++) {
    rc[i] = rc[i]/divisor;
  }
}



/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_typed_wts (int *type,
                        double *x,
                        double *y,
                        double *weight,
                        int *len,
                        int *typeA,
                        int *typeB,
                        double *r_low,
                        double *r,
                        int *len_r,
                        int *inds,
                        double *rc) {
  
  int i = 0;
  double divisor;
  double tmp_r_low = 0;
  double tmp_r = DBL_MAX;
  int tmp_len_r = 1;
  
  /*get the divisor in the pi function*/
  get_pi_typed_wts(type,x,y,weight,len,typeA,typeB,&tmp_r_low,&tmp_r,
                   &tmp_len_r,inds,&divisor);
  
  /*get the main pi function*/
  get_pi_typed_wts(type,x,y,weight,len,typeA,typeB,r_low,r,len_r,inds,rc);
  
  for (i = 0; i < *len_r; i++) {
    rc[i] = rc[i]/divisor;
  }
}



/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_clustsurvey (double *p,
                          double *x,
                          double *y,
                          int *s,
                          int *len,
                          double *r_low,
                          double *r,
                          int *len_r,
                          int *inds,
                          double *rc) {
  
  int i = 0;
  double divisor;
  double tmp_r_low = 0;
  double tmp_r = DBL_MAX;
  int tmp_len_r = 1;
  
  /*get the divisor in the pi function*/
  get_pi_clustsurvey(p,x,y,s,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
  
  /*get the main pi function*/
  get_pi_clustsurvey(p,x,y,s,len,r_low,r,len_r,inds,rc);
  
  for (i = 0; i < *len_r; i++) {
    rc[i] = rc[i]/divisor;
  }
}



/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_clustsurvey_wts (double *p,
                              double *x,
                              double *y,
                              int *s,
                              double *weight,
                              int *len,
                              double *r_low,
                              double *r,
                              int *len_r,
                              int *inds,
                              double *rc) {
  
  int i = 0;
  double divisor;
  double tmp_r_low = 0;
  double tmp_r = DBL_MAX;
  int tmp_len_r = 1;
  
  /*get the divisor in the pi function*/
  get_pi_clustsurvey_wts(p,x,y,s,weight,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
  
  /*get the main pi function*/
  get_pi_clustsurvey_wts(p,x,y,s,weight,len,r_low,r,len_r,inds,rc);
  
  for (i = 0; i < *len_r; i++) {
    rc[i] = rc[i]/divisor;
  }
}





/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_clustsurvey_wts2 (int *p,
                               double *x,
                               double *y,
                               double *s,
                               double *weight,
                               int *len,
                               double *r_low,
                               double *r,
                               int *len_r,
                               int *inds,
                               double *rc) {
  
  int i = 0;
  double divisor;
  double tmp_r_low = 0;
  double tmp_r = DBL_MAX;
  int tmp_len_r = 1;
  
  /*get the divisor in the pi function*/
  get_pi_clustsurvey_wts2(p,x,y,s,weight,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
  
  /*get the main pi function*/
  get_pi_clustsurvey_wts2(p,x,y,s,weight,len,r_low,r,len_r,inds,rc);
  
  for (i = 0; i < *len_r; i++) {
    rc[i] = rc[i]/divisor;
  }
}







