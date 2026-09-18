{..rs\Antillon\Dropbox\Harim\ESRFconSLS\ESRF_SLS_21_mayo_2020_20C.opa.txt}
{x=5.0E-03,px=0,y=1.0E-03,delta=0,100turns}

Energy = 3.000000;

BetaX   = 4.8032442; AlphaX  = 0.0000024;
EtaX    = 0.0027806; EtaXP   = 0.0000000;
BetaY   = 2.2903690; AlphaY  = 0.0000004;
EtaY    = 0.0000000; EtaYP   = 0.0000000;

{----- Variables ----------------------------------------------------}
ORBITDPP  = 0;
LSD       = 0.1;
F         = 0.8;
LSD1      = 0.1;
LSF1      = 0.1;
LSF1SH    = 0.05;
X1 =3.633167514008421e+00;
X2 =-4.258277621861492e+00;
X3 =-2.690860661253351e+00;
X4 =2.754505457254375e+00;
X5 =-3.336431720192718e+00;
X6 =-1.492176552197721e+00;
X7 =2.943728718874649e-03;
X8 =6.010287762639232e-01;
X9 =5.370545478298007e-01;
kse1=-31.67808026513335;
kfd2=-51.10107383195414;
kfd3=-93.324704323109;
ks1=-31.98668968024124;
ks2=41.64059482880812;
ksd3=-11.2882623813059;
ks1s=696.35209098583;
ks2s=488.0292055866639;
ko1=4.72312991538625;
ko2=46.6898786404116;
ko3=23.51531647572548;
ksf1=29.726909817;
ksd1=-127.01;
{----- FINV ----------------------------------------------------}
{----- Table of elements ----------------------------------------------------}
{----- se cambian estos elementos ----------------------------------------------------}
SE1    : Sextupole, L = LSD,        K =kse1,  N = 4, Ax = 50.00, Ay = 50.00;
FD2    : Sextupole, L = 0.094502,   K =kfd2,  N = 4, Ax = 50.00, Ay = 50.00;
FD3    : Sextupole, L = X7, K =kfd3,  N = 4, Ax = 50.00, Ay = 50.00;
S1     : Sextupole, L = LSD,        K =ks1,  N = 4, Ax = 50.00, Ay = 50.00;
S2     : Sextupole, L = LSD,        K =ks2,  N = 4, Ax = 50.00, Ay = 50.00;
SD3    : Sextupole, L = 0.010176,   K =ksd3,  N = 4, Ax = 50.00, Ay = 50.00;
S1S    : Sextupole, L = 0.002964,   K =ks1s,  N = 4, Ax = 50.00, Ay = 50.00;
S2S    : Sextupole, L = 0.172130,   K =ks2s,  N = 4, Ax = 50.00, Ay = 50.00;
O1     : Multipole, N = 4,          K =ko1, Ax = 50.00, Ay = 50.00;
O2     : Multipole, N = 4,          K =ko2, Ax = 50.00, Ay = 50.00;
O3     : Multipole, N = 4,          K =ko3, Ax = 50.00, Ay = 50.00;
SF1    : Sextupole, L = 2.2044e-01, K =ksf1,  N = 4, Ax = 50.00, Ay = 50.00;
SD1    : Sextupole, L = LSD,        K =ksd1,  N = 4, Ax = 50.00, Ay = 50.00;
{----- hasta aqui y la longitud de FD3   LFD3 ----------------------------------------------------}
{----- Table of segments ----------------------------------------------------}
D1     : Drift,      L = 2.654400-LSD, Ax = 5.00, Ay = 5.00;
D4     : Drift,      L = 0.081240, Ax = 5.00, Ay = 5.00;
D11    : Drift, L = 0.063628, Ax = 5.00, Ay = 5.00;
D12    : Drift, L = 0.0099526, Ax = 5.00, Ay = 5.00;
D5D6   : Drift, L =  X8, Ax = 5.00, Ay = 5.00;
D9D10  : Drift, L =  X9, Ax = 5.00, Ay = 5.00;
SF2    : Sextupole, L = 0.220440,   K =2.861070148054099e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
SD2    : Sextupole, L = LSD,        K =-6.959255479904094e+00*3,  N = 4, Ax = 50.00, Ay = 50.00;
SE2    : Sextupole, L = LSD,        K =1.218577437384940e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
S1SH   : Sextupole, L = 0.001379,   K =1.423515262859891e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
S3     : Sextupole, L = LSD,        K =-2.810527052262869e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
S4     : Sextupole, L = LSD,        K =-2.864350994439932e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
FD22   : Sextupole, L = 0.094502,   K =2.859707257217448e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
FD32   : Sextupole, L = 0.0017747,  K =1.317864852288966e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
SD32   : Sextupole, L = 0.010176,   K =1.140605742775726e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
S1SH2  : Sextupole, L = 0.001379,   K =1.083961208602719e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
S2S2   : Sextupole, L = 0.172130,   K =1.597721093016492e+01*3,  N = 4, Ax = 50.00, Ay = 50.00;
S1S2   : Sextupole, L = 0.002964,   K =6.184046275844920e+00*3,  N = 4, Ax = 50.00, Ay = 50.00;
O11    : Multipole, N = 4,          K =-4.661948609232911e+09*1e-8, Ax = 50.00, Ay = 50.00;
O21    : Multipole, N = 4,          K =-1.752515805260398e+09*1e-8, Ax = 50.00, Ay = 50.00;
O31    : Multipole, N = 4,          K =1.759193769861082e+10*1e-8, Ax = 50.00, Ay = 50.00;

DQ6    : Bending, L = 0.275390, T = -7.317925899999995E-01*F, K = 2.692600, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
A1     : Bending, L = 0.075497, T = 2.1719E-03*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
A2     : Bending, L = 0.384040, T = 5.3380E-01*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
A3     : Bending, L = 0.001995, T = 3.2534E-04*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
A4     : Bending, L = 0.913400, T = 2.0382E+00*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
A5     : Bending, L = 0.152490, T = 9.3133E-01*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B1     : Bending, L = 0.400570, T = 6.3294E-01*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B2     : Bending, L = 0.563170, T = 1.1254E+00*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B3     : Bending, L = 0.362720, T = 1.1741E+00*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B4     : Bending, L = 0.285610, T = 1.4465E+00*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B5     : Bending, L = 0.240960, T = 5.8358E-01*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B1S    : Bending, L = 0.015767, T = 8.0780E-02*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B2S    : Bending, L = 0.001644, T = -4.1155E-04*F, K = 0.000000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
B3S    : Bending, L = 0.212550, T = 1.7586E+00*F, K = 0.000000,  T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
DQ1S   : Bending, L = 0.257080, T = 8.1690E-01*F, K = -5.135300, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;
ABQ1   : Bending, L = 0.215990, T = -6.0542E-01*F, K = 6.191000, T1 = 0.000000, T2 = 0.000000, Ax = 50.00, Ay = 50.00;

QF1    : Quadrupole, L = 0.349140, K = X1, Ax = 5.00, Ay = 5.00;
QD2    : Quadrupole, L = 0.222950, K = X2, Ax = 5.00, Ay = 5.00;
QD3    : Quadrupole, L = 0.194780, K = X3, Ax = 5.00, Ay = 5.00;
QF4    : Quadrupole, L = 0.224580, K = X4, Ax = 5.00, Ay = 5.00;
QD5    : Quadrupole, L = 0.210950, K = X5, Ax = 5.00, Ay = 5.00;
QF7    : Quadrupole, L = 0.020986, K = X6, Ax = 5.00, Ay = 5.00;

DA1   : A1, A2, A3, A4, A5;
IDA1  : -DA1;
DB1 : B1, B2, B3, B4, B5;

DBA   : D1, SE1, QF1, FD2, QD2, FD3, IDA1, D4, QD3, SD1, O2, D5D6, S1, QF4, SF1,
O1, QF4,  S2, D9D10, O3, SD1, QD5, D11,  B1, B2, B3, B4, B5, D12, QF7, SD3, DQ6;

DBA2  : D1, SE2, QF1, FD22, QD2, FD32, IDA1, D4, QD3, SD2, O21, D5D6, S3, QF4, SF2, O11,
QF4, S4, D9D10, O31, SD2, QD5, D11,  B1, B2, B3, B4, B5, D12, QF7, SD32, DQ6;

IDBA  : -DBA;

CELA  :  S1S, ABQ1, S2S, DQ1S, B1S, B2S, B3S, B2S, B1S, DQ1S, S2S, ABQ1, S1S;
CELC  : S1SH, ABQ1, S2S, DQ1S, B1S, B2S, B3S, B2S, B1S, DQ1S, S2S2, ABQ1, S1SH2;
CELB  : S1SH2, ABQ1, S2S2, DQ1S, B1S, B2S, B3S, B2S, B1S, DQ1S, S2S2, ABQ1, S1S2;

CELD  : S1S, ABQ1, S2S, DQ1S, B1S, B2S, B3S, B2S, B1S, DQ1S, S2S, ABQ1, S1S;

CELL  : DBA, CELA, CELA, CELA, IDBA;
RING  : 20*CELL;

{..rs\Antillon\Dropbox\Harim\ESRFconSLS\ESRF_SLS_21_mayo_2020_20C.opa.txt}
