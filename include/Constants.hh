#ifndef BWCONSTANTS_H
#define BWCONSTANTS_H

struct Constants
{
    static constexpr double me   = 9.10938356e-31;
    static constexpr double eV   = 1.60217662e-19;
    static constexpr double c    = 2.99792458e8;
    static constexpr double hbar = 1.054571800e-34;
    static constexpr double kB   = 1.38064852e-23;
    static constexpr double pi   = 3.14159265359;
    static constexpr double r0   = 2.8179403227e-15;

    static constexpr double KeV   = 1.60217662e-16;
    static constexpr double MeV   = 1.60217662e-13;

    // Constants in MeV
    static constexpr double hbar_MeV = 6.582119514e-22;
    static constexpr double kb_MeV   = 8.6173303e-11;
    static constexpr double me_MeV   = 0.5109989461;
};

#endif
