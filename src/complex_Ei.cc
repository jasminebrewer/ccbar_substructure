/* routines for computing complex exponentials with high numerical precision are based on 
"On the Evaluation of the Complex-Valued Exponential Integral" by Vincent Pegoraro and Philipp Slusallek
 https://www.sci.utah.edu/~vpegorar/research/2011_JGT.pdf
 
 A test function is provided, which compares the computed values to results evaluated by Mathematica */

#include <complex>
#include <limits>
#include <iostream>
#include <cassert>
#include <iomanip>

using namespace std;

const std::complex<double> i(0, 1);
const double GAMMA = 0.57721566490153286060651209008; // euler-mascheroni constant
const int MAX_ITERATIONS = 100; // maximum iterations to try to obtain integral convergence
const double PRECISION = 1e-10; // value below which we consider a numerical quantity to be zero
const double MAX_ERROR = 1.e-5; // maximum error tolerance for integral convergence
const double DIST = 1.; // distance value determining the values of z for which we use the continued fraction representation. Value of 1 suggested by the above-mentioned paper
const double INF = numeric_limits<float>::infinity();

double sign(double x){
    if (x > 0.) return 1.;
    if (x < 0.) return -1.;
    return 0.;
}

bool converged(complex<double> ei, complex<double> old) {
    
    return abs((ei - old).real()) <= MAX_ERROR * abs(ei.real()) && abs((ei - old).imag()) <= MAX_ERROR * abs(ei.imag());
}

// compute Ei using the Asymptotic Series
std::complex<double> EiAsymptoticSeries(complex<double> z) {

    complex<double> ei = sign(z.imag()) * i * M_PI;
    complex<double> tmp = exp(z) / z;
    complex<double> old;
    for (int k=1; k<=floor(abs(z))+1; k++) {
        old = ei;
        ei += tmp;
        if (converged(ei, old)) return ei;
        tmp *= double(k) / z;
    }
    return INF;
}

// compute Ei using the Power Series
std::complex<double> EiPowerSeries(std::complex<double> z) {

    std::complex<double> ei = GAMMA + log(abs(z)) + sign(z.imag()) * i * abs(arg(z));
    std::complex<double> tmp = 1.;
    std::complex<double> old;
    for (int k=1; k<=MAX_ITERATIONS; k++) {
        tmp *= z / double(k);
        old = ei;
        ei += tmp / double(k);
        if (converged(ei, old)) return ei;
    }
    return INF;
}


// compute Ei using the Continued Fraction representation
std::complex<double> EiContinuedFraction(std::complex<double> z) {

    std::complex<double> ei = sign(z.imag()) * i * M_PI;
    std::complex<double> c, d, old;
    if ( abs(ei) > PRECISION ) { // ei is not zero
        c = 1. / ei;
        d = 0.;
        c = 1. / (1. - z - exp(z) * c);
        d = 1. / (1. - z - exp(z) * d);
        ei *= d / c;
    }
    else {
        c = INF;
        d = 0.;
        c = 0.;
        d = 1. / (1. - z - exp(z) * d);
        ei = d * (- exp(z));
    }
    for(int k=1; k<=MAX_ITERATIONS; k++) {
        c = 1. / (2. * double(k) + 1. - z - double(k) * double(k) * c);
        d = 1. / (2 * double(k) + 1. - z - double(k) * double(k) * d);
        old = ei;
        ei *= d / c;
        if (converged(ei, old)) return ei;
    }
    return INF;
}

// primary function for evaluating Ei, that decides which representation to use to compute Ei for a given value of z.
std::complex<double> Ei(std::complex<double> z) {

    std::complex<double> ei;
    if (abs(z) > 1e10) return sign(z.imag()) * i * M_PI + exp(z) / z;
    if (abs(z) > 2 - 1.035 * log(MAX_ERROR)) {
        ei = EiAsymptoticSeries(z);
        if (ei != INF) return ei;
    }
    if (abs(z) > DIST && (z.real() < 0 || abs(z.imag()) > DIST)) {
        ei = EiContinuedFraction(z);
        if (ei != INF) return ei;
    }
    if (abs(z) > 0) {
        ei = EiPowerSeries(z);
        if (ei != INF) return ei;
    }
    if (abs(z) < PRECISION) return -INF;
}

// testing function for the Ei function given a wide range of inputs z=a+bi, with a,b \in (0,10).
// numerical outputs are compared to results from Mathematica
void test_Ei() {

    // the percent (relative to the absolute value) of numerical tolerance accepted for the agreement of the solutions
    double tolerance = 1.e-5;
    
    std::complex<double> input = 0.+0.1*i;
    std::complex<double> answer = -1.727868387 + 1.670740788*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );
 
    input = 0.+0.5*i;
    answer = -0.1777840788 + 2.063903745*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.+1.*i;
    answer = 0.3374039229 + 2.516879397*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.+2.*i;
    answer = 0.4229808288 + 3.176209304*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.+5.*i;
    answer = -0.1900297497 + 3.120727572*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.+10.*i;
    answer = -0.045456433 + 3.229143921*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+0.*i;
    answer = -1.622812814;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+0.1*i;
    answer = -1.278911182 + 0.890509206*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+0.5*i;
    answer = -0.05989028577 + 1.891825029*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+1.*i;
    answer = 0.4284756357 + 2.464724315*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+2.*i;
    answer = 0.4702003722 + 3.19928769*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+5.*i;
    answer = -0.2101365348 + 3.114556763*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.1+10.*i;
    answer = -0.05122261343 + 3.237939129*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+0.*i;
    answer = 0.4542199049;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+0.1*i;
    answer = 0.4703190328 + 0.3270586372*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+0.5*i;
    answer = 0.7139426606 + 1.424048023*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+1.*i;
    answer = 0.9228991906 + 2.32587992*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+2.*i;
    answer = 0.6935676688 + 3.346793592*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+5.*i;
    answer = -0.3119980112 + 3.077360022*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 0.5+10.*i;
    answer = -0.08218694055 + 3.282571009*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+0.*i;
    answer = 1.895117816;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+0.1*i;
    answer = 1.895095329 + 0.2713771604*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+0.5*i;
    answer = 1.883148319 + 1.307944759*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+1.*i;
    answer = 1.764625986 + 2.387769852*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+2.*i;
    answer = 1.042167708 + 3.701501426*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+5.*i;
    answer = -0.5031006696 + 2.987318243*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 1.+10.*i;
    answer = -0.1468901084 + 3.367311464*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 2.+0.*i;
    answer = 4.954234356;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 2.+2.*i;
    answer = 1.892078162 + 5.316962438*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 5.+2.*i;
    answer = 3.230031092 + 36.30913195*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 5.+10.*i;
    answer = -11.13434095 + 11.08683318*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 10.+0.5*i;
    answer = 2248.725238 + 1064.093276*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 10.+2.*i;
    answer = -485.0505411 + 2374.818989*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    input = 10.+10.*i;
    answer = -1576.150427 + 436.9192317*i;
    assert ( (abs(Ei(input) - answer) < tolerance*abs(answer)) || printf("failed assert for  %f%+fi: your answer is %f%+fi while the test is %f%+fi \n", input.real(), input.imag(), Ei(input).real(), Ei(input).imag(), answer.real(), answer.imag()) );

    std::complex<double> x = 10.+10.*i;
    std::complex<double> y = 1.001 + 0.001*i;
    cout << "diff: " << Ei(x)-Ei(x/y) << endl;

    cout << "completed Ei unit tests successfully, with all errors satisfying abs(error) < " << tolerance << "*abs(answer)" << endl;

}