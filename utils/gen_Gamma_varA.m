function Gamma = gen_Gamma_varA(alpha,P)
    B_m = TriaToepMulShort(alpha,P,false);
    C_m = TriaToepMulShort(alpha,P,true);
    Gamma = 1/alpha(1) * ((B_m - C_m));
end