function [rIIII,E,resG, resL] = OverallRateOfChange_rxn(pA, CB)

global b c DiffA DiffB HA KLa KGa
     
% codice completo

    %% Preliminary calculations
    resG = 1/KGa;
    
    %% Enhancement factor
    % Region III

    pAI = pA;
        
    for i=1:100
        Ei = 1+ DiffB/DiffA*CB*HA/b/pAI;
        E =Ei;

        resL = HA/KLa/E;
        rIIII = pA/(resG+resL); 
        rIIIIGas = KGa*(pA-pAI);
        if(abs(rIIII-rIIIIGas)/rIIII<0.001)
            break;
        end
        pAI = pA-rIIII/KGa;  
    end  
    
end
