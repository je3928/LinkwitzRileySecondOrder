//
//  LinkwitzRiley2ndOrder.h
// 
//  Author: Jordan Evans
//  22.03.2022
//  2nd Order Linkwitz-Riley filters for use in cross over networks. Adapted from earlevel engineerings biquad object.

// Update 15.12.2022: filter structure changed to direct form 2 transposed 

// Notes: type 0 = LPF, type 1 = HPF

#include <math.h>



class LinkwitzRiley2ndOrder {
public:
    LinkwitzRiley2ndOrder();
    LinkwitzRiley2ndOrder(int type, double Fc, double Fs);
    ~LinkwitzRiley2ndOrder();
    void setType(int type);
    void setFc(double Fc);
    void setFilter(int type, double Fc, double Fs);
    float process(float in);

    double pi = 3.14159265358979323846;

protected:
    void calcFilter(void);

    int type;
    double a0, a1, a2, b1, b2;
    double Fc, Fs;
    double s1, s2, s1minus1, s2minus1;
};

inline float LinkwitzRiley2ndOrder::process(float in) {
    double out = a0 * in + s1minus1;
    s1 = s2minus1 + a1 * in - b1 * out;
    s2 = a2 * in - b2 * out;


    s1minus1 = s1;
    s2minus1 = s2;
    return out;
}
