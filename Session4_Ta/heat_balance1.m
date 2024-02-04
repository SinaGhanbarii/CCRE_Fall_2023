function F=heat_balance1(y)
    global V c0 c_0 Q rho cp delH_rxn T_feed nu UA k0 Tcool
    T=y(1);
    c=[y(2) y(3)];
    
    k=k0*exp(-12000/T);
    %Part 1
    F(1)=((Q*c0*delH_rxn)/( Q/(k*V)+1))-(Q*rho*cp*(T-T_feed))-UA*(T-Tcool)-5000*(T-280);
    %Part 2
    %F(1)=((Q*c0*delH_rxn)/( Q/(k*V)+1))-(Q*rho*cp*(T-T_feed))-UA*(T-Tcool);%-5000*(T-280);

    c(1)=max(c(1),1e-6);
    
    
    %reaction rate
    
    r(1)=k*c(1);
    
    R=zeros(2,1);
    
    
    for i=1:2
        R(i)=0;
        R(i)=R(i)+nu(i)*r;
    end
    
    for i=1:2
        F(i+1)=Q*c_0(i)-Q*c(i)+R(i)*V;
    end
    end
    
    
    
    
    
    
    
    
    
    
    
    