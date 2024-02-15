function [rIIII,E, Ei, MH, resG, resL, resK] = OverallRateOfChange_KN(pA, CB)

global b c DiffA DiffB fL HA a KLa KGa kl 
     
% codice completo

    %% Preliminary calculations
    klStar = kl*CB;
    klPrime = klStar*CB;
    KL = KLa/a;
    MH = sqrt(klPrime*DiffA)/KL;
    resG = 1/KGa;
    resK = HA/(klPrime*fL);
    
    %% Enhancement factor
    if (MH<1)               % Regions I and II
        
      if (MH<0.3)
          E = 1;
      else 
          E = 1+MH^2/3;
      end

      resL = HA/KLa/E;
      rIIII = pA/(resG+resL+resK);

    else                    % Region III

        pAI = pA;
        
        for i=1:100
            Ei = 1+ DiffB/DiffA*CB*HA/b/pAI;

            if (MH<5*Ei)
                E = MH;
            else
                E =Ei;
            end
            
            resL = HA/KLa/E;
            rIIII = pA/(resG+resL+resK);
            
            rIIIIGas = KGa*(pA-pAI);
            
            if(abs(rIIII-rIIIIGas)/rIIII<0.001)
                break;
            end
            
            pAI = pA-rIIII/KGa;  
        end  
    end
    
end
