//
//  LinkwitzRiley2ndOrder.h
// 
//  Author: Jordan Evans
//  22.03.2022
//  2nd Order Linkwitz-Riley filters for use in cross over networks. Adapted from earlevel engineerings biquad object.

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
    double z1, z2;
};

inline float LinkwitzRiley2ndOrder::process(float in) {
    double out = in * a0 + z1;
    z1 = in * a1 + z2 - b1 * out;
    z2 = in * a2 - b2 * out;
    return out;
}
