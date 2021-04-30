/*Random number generators. ranf0 is to be used any time we need a uniform random
number, excluding as a test number for a Gaussian generator.  ranf1,2,3 are to be 
used with gaussian1,2,3, respectively.*/


double ranf0()
/* Uniform random number generator x(n+1)= a*x(n) mod c
   with a = pow(7,5) and c = pow(2,31)-1.
   Copyright (c) Tao Pang 1997. */
{
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
int il,ih,it;
double rc;
extern int iseed0;
ih = iseed0/iq;
il = iseed0%iq;
it = ia*il-ir*ih;
if (it > 0)
  {
  iseed0 = it;
  }
else
  {
iseed0 = ic+it;
  }
rc = ic;

return iseed0/rc;
}//end of ranf0()


double ranf2()
/* Uniform random number generator x(n+1)= a*x(n) mod c
   with a = pow(7,5) and c = pow(2,31)-1.
   Copyright (c) Tao Pang 1997. */
{
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
int il,ih,it;
double rc;
extern int iseed2;
ih = iseed2/iq;
il = iseed2%iq;
it = ia*il-ir*ih;
if (it > 0)
  {
  iseed2 = it;
  }
else
  {
iseed2 = ic+it;
  }
rc = ic;
return iseed2/rc;
}//end of ranf2()



double ranf4()
/* Uniform random number generator x(n+1)= a*x(n) mod c
   with a = pow(7,5) and c = pow(2,31)-1.
   Copyright (c) Tao Pang 1997. */
{
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
int il,ih,it;
double rc;
extern int iseed4;
ih = iseed4/iq;
il = iseed4%iq;
it = ia*il-ir*ih;
if (it > 0)
  {
  iseed4 = it;
  }
else
  {
iseed4 = ic+it;
  }
rc = ic;
return iseed4/rc;
}//end of ranf4()


double ranf6()
/* Uniform random number generator x(n+1)= a*x(n) mod c
 *    with a = pow(7,5) and c = pow(2,31)-1.
 *       Copyright (c) Tao Pang 1997. */
{
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
int il,ih,it;
double rc;
extern int iseed6;
ih = iseed6/iq;
il = iseed6%iq;
it = ia*il-ir*ih;
if (it > 0)
  {
  iseed6 = it;
  }
else
  {
iseed6 = ic+it;
  }
rc = ic;
return iseed6/rc;
}//end of ranf6()




/*Takes uniformly distributed random numbers and gives random number according to
Gaussian distribution by the rejection method.*/

//rotational brownian motion
double gaussian2()
{
       double y, py, ptest;      
 
       y = gauss_range2*(ranf2() - 0.5);
       ptest = gauss_prefact2*ranf2();//gauss_prefact is the maximum value of P_y(y), the probability density function
       
       py = gauss_prefact2*exp(gauss_exp_const2*y*y);
       
       if(ptest < py)
          return(y);
       else
          return gaussian2();
}




///////////////////////////////////////////////////////////////////
double ranf5()
/* Uniform random number generator x(n+1)= a*x(n) mod c
 *    with a = pow(7,5) and c = pow(2,31)-1.
 *       Copyright (c) Tao Pang 1997. */
{
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
int il,ih,it;
double rc;
extern int iseed5;
ih = iseed5/iq;
il = iseed5%iq;
it = ia*il-ir*ih;
if (it > 0)
  {
  iseed5 = it;
  }
else
  {
iseed5 = ic+it;
  }
rc = ic;

return iseed5/rc;
}//end of ranf5()


//gaussian random number generator for heavy polymer mono rotation
double gaussian6()
{
       double y, py, ptest;

       y = gauss_range6*(ranf6() - 0.5);
       ptest = gauss_prefact6*ranf6();//gauss_prefact is the maximum value of P_y(y), the probability density function

       py = gauss_prefact6*exp(gauss_exp_const6*y*y);

       if(ptest < py)
          return(y);
       else
          return gaussian6();
}


/*****************/
double gaussian_std()
{
       double y, py, ptest;

       y = gauss_range_std*(ranf_std() - 0.5);
       ptest = ranf_std();//gauss_prefact is the maximum value of P_y(y), the probability density function
       py = exp(gauss_exp_const_std*y*y);


       if(ptest < py)
          return(y);
       else
          return gaussian_std();
}



double gaussian_inverf(unsigned long long step)
{
	double r;

	do{
	  r = RNUM0.get_double();
	}while(r < 1.e-16);

	return (SQRT_TWO* erf_inv( 2.*(r-0.5) ));
}

