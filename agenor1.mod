//variables under investigation
var cki chaf cxf cgrowth cefr cefw cefp crtotal;

parameters cwi cvi ctax cb cbeta cmui cthetar cpc cxr cv3 ck2 cvp calpha cwh cvh cmuh cva ck3 cwe cve cv1 cv2 cmconst cnq cpiq czetap cnn cxbar cgammab cpa crho cnfh cnmh cnfc cnmc cmub cxx cnh cnc csig cnt ck1;

//exogenous calibrated parameter values
calpha = 0.15;
cwh = 0.39;
cvh = 0.115;
cmuh = 0.8;
cva = 0.2;
ck3 = 0.2;
cwe = 0.39;
cve = 0.171;
cv1 = 0.4;
cv2 = 0.15;
cmconst = 2.011;
cnq = 7.634;
cpiq = 0.7;
czetap = 1.0;
cnn = 4.289;
cxbar = 0.123;
cwi = 0.39;
cvi = 0.05;
ctax = 0.239;
cb = 0.6;
cbeta = 0.35;
cmui = 1;
cthetar = 0.105;
cpc = 0.854;
cxr = 0.6;
cv3 = 0.1;
ck2 = 0.15;
cvp = 0.8;
cgammab = 0.5;
cpa = 0.702;
crho = 0.04;
cnmh = 6.0;
cnmc = 4.5;
cnfc = 1.845;
cnfh = 10.142;
cmub = 1.0;
ck1 = 0.5;
//cperiod = number of years in a period
cperiod = 20;

//endogenously calibrated parameters
cxx = (cxbar^(1 - cgammab)) * (cxr / (1 - cxr))^(-1 * cv3 * cmub * cgammab);
cnc = cxx * cnfc + (1 - cxx) * cnmc;
cnh = cxx * cnfh + (1 - cxx) * cnmh;
csig = cpa / ((1 + crho)^cperiod * cnc + cpa);
cnt = (1 - (cnh * ck2) / cnn) / ((1 - ((cnh * ck2) / cnn) + (cnc / (cnn * (1 - csig)))) * cthetar * cpc);


model;

//public-private capital ratio
cki = (((cwi * cvi * ctax * (1 + cb) * cbeta)^(cmui) / (cb * cbeta * (1 - ctax) * (cb^(-1) + 1) * csig * (1 - cthetar * cpc * cnt) * ((cxr / (1 - cxr))^(cbeta * (cv3 + ck2 * cvp)))^(1 - cmui))) * (cefw)^(-cbeta * (1 - cmui)) * ((chaf)^(cvp) / cxf)^(-2 * cbeta * (1- cmui)))^(1 / (1 - (1 - cmui) * (1 - calpha)));

//female health status
chaf = ((1 - cxr)^ck2 * ((cwh * cvh * ctax * (1 + cb) * cbeta)^cmuh * ((cxr / (1 - cxr))^(cbeta * (cv3 + ck2 * cvp)))^cmuh)^ck3 * (cxr / (1 - cxr))^(-cv3 * cva) * cefr^ck2 * cki^(ck3 * (1 - cmuh * (1 - calpha))) * cefw^(ck3 * cbeta * cmuh) * cxf^(-2 * ck3 * cbeta * cmuh))^(1 / (1 - ck1 - cvp * 2 * ck3 * cbeta * cmuh));

//private capital effective labour j ratio
cxf = ((cb * cbeta * (1 - ctax) * (cb^(-1) + 1) * csig * (1 - cthetar * cpc * cnt)) / ((1 - cxr)^cv3 * (cpc * cnt)^(1 - cv1) * 0.5) * (cwe * cve * ctax * (1 + cb) * cbeta)^(-cv1) * ((cxr / (1 - cxr))^(cbeta * (cv3 + ck2 * cvp)))^(1 - cv1) * cki^(-cv2 + calpha * (1 - cv1)) * chaf^(2 * cvp * cbeta * (1 - cv1)) * cefr^(-cv3) * cefw^(cbeta * (1 - cv1)))^(1 / (1 - (1 - 2 * cbeta) * (1 - cv1)));

//growth rate of final output
cgrowth = cmconst * (cxr / (1 - cxr))^(cbeta * (cv3 + ck2 * cvp)) * cki^calpha * cefw^cbeta *(cbeta * csig * (1 - cthetar * cpc * cnt)) / (((1 - ctax) * (1 + cb))^(-1)) * (chaf)^(cvp * 2 * cbeta) * cxf^(-2 * cbeta);

//time spent on home production
cefp = (1 + ((cnq * cpiq * (1 - csig) * cnc^(-1)) / (1 + cnh * ck2 * (1 - csig) * cnc^(-1))))^-1 * ((cnq * cpiq * (1 - csig) * cnc^(-1)) / (1 + cnh * ck2 * (1 - csig) * cnc^(-1)) - czetap * cki);

//time allocated to each child
cefr = (((1 - ((cnh * ck2) / cnn) + (cnc / (cnn * (1 - csig)))) * cthetar * cnh * ck2 * (1 - csig)) / (cnc * (1 - (cnh * ck2 / cnn)))) * ((1 - cefp) / (1 + cnh * ck2 * (1 - csig) * cnc^(-1)));

//time spent on market work
cefw = 1 - cefp - cpc * cnt * cefr;

//time allocated to all children
crtotal = cpc * cnt * cefr;

end;


steady;





