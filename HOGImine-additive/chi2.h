#ifndef CCHI2
#define CCHI2

#include<math.h>


double regularizedLowerIncompleteGamma(double x, double alpha) {
	if(x <= 0.0 || alpha <= 0.0) return 0.0;

    double gamma_f = exp(alpha*log(x) - lgamma(alpha) - x);
	if(x < alpha + 1.0) {
		double i = alpha;
		double tmp_sum = 1.0;
		double sum = tmp_sum;
		while(tmp_sum/sum > 1e-10) {
			i++;
			tmp_sum *= x/i;
			sum += tmp_sum;
		}
		return gamma_f*sum/alpha;
	} else {
        //Solve via evaluting continued fractions
		//Evaluation of continued fraction
		double a=1.0-alpha;
		double b=1+x+a;
		double pa1 = 1.0;
		double pb1 = x;
		double pa2 = x + 1.0;
		double pb2 = b*x;
		double func = pa2/pb2;
		double pa,pb,ratio,tmp;
		double i = 0;

		while(1) {
			i++;
			a++;
			b += 2.0;
			pa = b*pa2-a*i*pa1;
			pb = b*pb2-a*i*pb1;
			if(pb) {
				ratio = pa/pb;
				tmp = fabs((func-ratio));
				if(tmp<=1e-10*ratio) break;
				func = ratio;
			} else tmp=1.0;
			pa1=pa2;
			pb1=pb2;
			pa2=pa;
			pb2=pb;
			if(i>100) break;//Maximum number if iterations
		}
		return 1.0-func*gamma_f;
	}
}

/*
Function to compute a complemented incomplete gamma function based on a continued fraction
*/
double complementedIncompleteGamma(double x, double alpha) {
	if((x <= 0) || ( alpha <= 0))
		return 1.0;

    if((x < 1.0) || (x < alpha))
		return 1.0 - regularizedLowerIncompleteGamma(x,alpha);

    double gamma_f = exp(alpha*log(x) - lgamma(alpha) - x);

	//continued fraction
	double y = 1.0 - alpha;
	double z = 1.0 + x + y;
	double c = 0.0;
	double pkm2 = 1.0;
	double qkm2 = x;
	double pkm1 = x + 1.0;
	double qkm1 = z * x;
	double func = pkm1/qkm1;
	double i = 0;
	double ratio,tmp,pk,qk;

	while(1) {
		i++;
		c += 1.0;
		y += 1.0;
		z += 2.0;
		pk = pkm1 * z  -  pkm2 * y*c;
		qk = qkm1 * z  -  qkm2 * y*c;
		if( qk != 0 ) {
			ratio = pk/qk;
			tmp = fabs( (func - ratio)/ratio );
			if(tmp<=1e-10*ratio) break;
			func = ratio;
		} else {
			tmp = 1.0;
		}
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;
		if( fabs(pk) > 1e32) {
			pkm2 *= 1e-32;
			pkm1 *= 1e-32;
			qkm2 *= 1e-32;
			qkm1 *= 1e-32;
		}
		if(i>100) break; //Max number of iterations
	}
	return func * gamma_f;
}

/*
Computed the survival function of the chi2 distribution function for x and k
*/
inline double Chi2_sf(double x, double k) {
    return complementedIncompleteGamma(0.5*x,0.5*k);
}

/*
Computed the cumulative distribution function of the chi2 distribution function for x and k
*/
inline double Chi2_cdf(double x, double k) {
    if (k==2.0) {
        return 1.0 - exp(-0.5*x);
    } else {
        return regularizedLowerIncompleteGamma(0.5*x,0.5*k);
    }
}

/*
Probability density function of Chi2
*/
inline double Chi2_pdf(double x, double k) {
    if (x<0.0) return 0.0;
    return pow(x,0.5*k-1.0)*exp(-0.5*x)/(pow(2.0,0.5*k)*tgamma(0.5*k));
}

#endif
