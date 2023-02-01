%% Set up the right-hand side of the system (1) in the manuscript

function dy=odeSystem_singleRho(~,y,p)

   n=y(1); a=y(2); m1=y(3); m2=y(4); c=y(5); g=y(6);
   h=y(7); ha=y(8); e=y(9); s=y(10); sa=y(11); 
   
   dy(1)= p.rho*c/(1+g)-p.nu*(1+g/p.betaG)*n/(1+c/p.betaC);
   dy(2)= p.nu*(1+g/p.betaG)*n/(1+c/p.betaC)-p.rho*p.gammaA*a-p.phi*a*(m1+p.phi2*m2);
   dy(3)= p.rho*c-p.rho*p.km1*p.phi*a*m1+p.km2*m2-p.gammaM*m1*(1+p.gammaM2*m2);
   dy(4)= p.rho*p.km1*p.phi*a*m1-p.km2*m2-p.gammaM*m2;
   dy(5)= p.rho*p.kn*n^2/(p.betaN^2+n^2)+ p.rho*p.gammaA*(a^2./(p.betaA^2+a^2))+ha+p.km*m1-c;
   dy(6)= p.rho*p.kg*m2+p.kh*h-p.gammaG*g;
   dy(7)= p.chiH*p.phi*s*ha*(m1+p.phi2*m2)+p.gammaE*e-p.nu2*h*c+p.gammaH*ha*s;
   dy(8)= p.nu2*h*c-p.chiH*p.phi*ha*(m1+p.phi2*m2)-p.gammaH*ha;
   dy(9)= p.chiH*p.phi*sa*ha*(m1+p.phi2*m2)-p.gammaE*e+p.gammaH*ha*sa;
   dy(10)= p.r2*sa*(1+g)-p.rho*p.r1*s*c;
   dy(11)= p.rho*p.r1*s*c-p.r2*sa*(1+g);
   
   dy=dy';

end
