%Apply Horner's algorithm to calculate (probabilist's) Hermite polynomials
function [res] = probHermiteH(n,x)
if n<=5
    switch n
    case 0
        res = 1;
    case 1
        res = x;
    case 2
        res = x*x - 1;
    case 3
        res = x * (x*x - 3);
    case 4
        x2 = x*x;
        res = 3 + x2 * (-6 + x2);
    case 5
        x2 = x*x;
        res = x * (15 + x2 * (-10 + x2));
    case 6
        x2 = x*x;
        res = -15 + x2 * (45 + x2 * (-15 + x2));
    case 7
        x2 = x*x;
        res = x *(-105 + x2 * (105 + x2 *(-21 + x2)));
    case 8
        x2 = x*x;
        res = 105 + x2 *(-420 + x2 *(210 + x2 *(-28 + x2)));
    case 9
        x2 = x*x;
        res = x *(945 + x2 *(-1260 + x2 *(378 + x2 *(-36 + x2))));
    case 10
        x2 = x*x;
        res = -945 + x2 *(4725 + x2 *(-3150 + x2 *(630 + x2 *(-45 + x2))));
    case 11
        x2 = x*x;
        res = x *(-10395 + x2 *(17325 + x2 *(-6930 + x2 *(990 + x2 *(-55 + x2)))));
    case 12
        x2 = x*x;
        res = 10395 + x2 *(-62370 + x2 *(51975 + x2 *(-13860 + x2 *(1485 + x2 *(-66 + x2)))));
    case 13
        x2 = x*x;
        res = x *(135135 + x2 *(-270270 + x2 *(135135 + x2 *(-25740 + x2 *(2145 + x2 *(-78 + x2))))));
    case 14
        x2 = x*x;
        res = -135135 + x2 *(945945 + x2 *(-945945 + x2 *(315315 + x2 *(-45045 + x2 *(3003 + x2 *(-91 + x2))))));
    case 15
        x2 = x*x;
        res = x *(-2027025 + x2 *(4729725 + x2 *(-2837835 + x2 *(675675 + x2 *(-75075 + x2 *(4095 + x2 *(-105 + x2)))))));
    end
else
    res = x * probHermiteH(n-1,x) - (n-1) * probHermiteH(n-2,x);
end