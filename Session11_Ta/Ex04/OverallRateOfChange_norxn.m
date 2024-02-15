function rIIII = OverallRateOfChange_norxn(pA, cA)

global HA KLa KGa 
     
% codice completo

    %% Preliminary calculations
    resG = 1/KGa;
    resL = HA/(KLa);
    rIIII = (pA-HA*cA)/(resG+resL);
end
