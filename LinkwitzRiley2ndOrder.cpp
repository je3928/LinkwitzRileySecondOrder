//
//  LinkwitzRiley2ndOrder.cpp
// 
//  Author: Jordan Evans
//  22.03.2022
//  2nd Order Linkwitz-Riley filters for use in cross over networks. Adapted from earlevel engineerings biquad object.

// Update 15.12.2022: filter structure changed to direct form 2 transposed 

// Notes: type 0 = LPF, type 1 = HPF



#include "LinkwitzRiley2ndOrder.h"

LinkwitzRiley2ndOrder::LinkwitzRiley2ndOrder() {
    type = 0;
    a0 = 1.0;
    a1 = a2 = b1 = b2 = 0.0;
    Fc = 0.50;
    Fs = 0.0000000;
    s1 = s2 = s1minus1 = s2minus1 = 0.0;
}

LinkwitzRiley2ndOrder::LinkwitzRiley2ndOrder(int type, double Fc, double Q) {
    setFilter(type, Fc, Q);
    s1 = s2 = s1minus1 = s2minus1 = 0.0;
}

LinkwitzRiley2ndOrder::~LinkwitzRiley2ndOrder() {
}

void LinkwitzRiley2ndOrder::setType(int type) {
    this->type = type;
    calcFilter();
}


void LinkwitzRiley2ndOrder::setFc(double Fc) {
    this->Fc = Fc;
    calcFilter();
}

void LinkwitzRiley2ndOrder::setFilter(int type, double Fc, double Fs) {

    s1 = s2 = s1minus1 = s2minus1 = 0.0;

    this->type = type;
    this->Fs = Fs;
    this->Fc = Fc;
    calcFilter();
}

void LinkwitzRiley2ndOrder::calcFilter(void) {

    double theta = pi * Fc / Fs;
    double omega = pi * Fc;
    double kappa = omega / tan(theta); 
    double delta = (kappa * kappa) + (omega * omega) + (2 * kappa * omega);
    switch (this->type) {
    case 0:
        a0 = (omega*omega) / delta;
        a1 = 2 * ((omega*omega) / delta);
        a2 = (omega*omega) / delta;
        b1 = ((-2 * (kappa*kappa)) + (2*(omega * omega))) / delta;
        b2 = ((-2 * kappa * omega) + (kappa * kappa) + (omega * omega)) / delta;
        break;

    case 1:
        a0 = (kappa * kappa) / delta;
        a1 = (-2 * (kappa * kappa)) / delta;
        a2 = (kappa * kappa)/delta;
        b1 = ((-2 * (kappa * kappa)) + (2*(omega * omega))) / delta;
        b2 = ((-2 * kappa * omega) + (kappa * kappa) +  (omega * omega)) / delta;
        break;
    }

    return;
}