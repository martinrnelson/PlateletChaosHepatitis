function dy=odeSystemLyapunov_manyRhos(~,y,p)

   n=y(1); a=y(2); m1=y(3); m2=y(4); c=y(5); g=y(6);
   h=y(7); ha=y(8); s=y(9);
   
   % Equations
   dy(1)= p.rho(1)*c/(1+g)-p.nu*(1+g/p.betaG)*n/(1+c/p.betaC);
   dy(2)= p.nu*(1+g/p.betaG)*n/(1+c/p.betaC)-p.rho(2)*p.gammaA*a-p.phi*a*(m1+p.phi2*m2);
   dy(3)= p.rho(3)*c-p.rho(4)*p.km1*p.phi*a*m1+p.km2*m2-p.gammaM*m1*(1+p.gammaM2*m2);
   dy(4)= p.rho(4)*p.km1*p.phi*a*m1-p.km2*m2-p.gammaM*m2;
   dy(5)= p.rho(5)*p.kn*n^2/(p.betaN^2+n^2)+ p.rho(2)*p.gammaA*(a^2./(p.betaA^2+a^2))+ha+p.km*m1-c;
   dy(6)= p.rho(6)*p.kg*m2+p.kh*h-p.gammaG*g;
   dy(7)= p.chiH*p.phi*s*ha*(m1+p.phi2*m2)+p.gammaE*(1-h-ha)-p.nu2*h*c+p.gammaH*ha*s;
   dy(8)= p.nu2*h*c-p.chiH*p.phi*ha*(m1+p.phi2*m2)-p.gammaH*ha;
   dy(9)= p.r2*(1-s)*(1+g)-p.rho(7)*p.r1*s*c;
   
   Y=reshape(y(10:end),9,9);
   
   % Jacobian
   Jac=zeros(9);
   Jac(1,1)=-(p.nu*(g/p.betaG + 1))/(c/p.betaC + 1);
   Jac(1,5)=p.rho(1)/(g + 1) + (n*p.nu*(g/p.betaG + 1))/(p.betaC*(c/p.betaC + 1)^2);
   Jac(1,6)=- p.rho(1)*c/(g + 1)^2 - (n*p.nu)/(p.betaG*(c/p.betaC + 1));
   Jac(2,1)=(p.nu*(g/p.betaG + 1))/(c/p.betaC + 1);
   Jac(2,2)=-p.gammaA*p.rho(2) - p.phi*(m1 + m2*p.phi2);
   Jac(2,3)=-p.phi*a;
   Jac(2,4)=-p.phi*p.phi2*a;
   Jac(2,5)=-(n*p.nu*(g/p.betaG + 1))/(p.betaC*(c/p.betaC + 1)^2);
   Jac(2,6)=(n*p.nu)/(p.betaG*(c/p.betaC + 1));
   Jac(3,2)=-p.rho(4)*p.km1*p.phi*m1;
   Jac(3,3)=-p.gammaM*(p.gammaM2*m2 + 1) - a*p.rho(4)*p.km1*p.phi;
   Jac(3,4)=p.km2 - p.gammaM*p.gammaM2*m1;
   Jac(3,5)=p.rho(3);
   Jac(4,2)=p.rho(4)*p.km1*p.phi*m1;
   Jac(4,3)=a*p.rho(4)*p.km1*p.phi;
   Jac(4,4)=-p.gammaM-p.km2;
   Jac(5,1)=p.rho(5)*p.kn*((2*n)/(p.betaN^2 + n^2) - (2*n^3)/(p.betaN^2 + n^2)^2);
   Jac(5,2)=p.rho(2)*p.gammaA*((2*a)/(a^2 + p.betaA^2) - (2*a^3)/(a^2 + p.betaA^2)^2);
   Jac(5,3)=p.km;
   Jac(5,5)=-1;
   Jac(5,8)=1;
   Jac(6,4)=p.rho(6)*p.kg;
   Jac(6,6)=-p.gammaG;
   Jac(6,7)=p.kh;
   Jac(7,3)=p.chiH*ha*p.phi*s;
   Jac(7,4)=p.chiH*ha*p.phi*p.phi2*s;
   Jac(7,5)=-p.nu2*h;
   Jac(7,7)=-p.gammaE-p.nu2*c;
   Jac(7,8)=p.gammaH*s - p.gammaE + p.chiH*p.phi*s*(m1 + m2*p.phi2);
   Jac(7,9)=p.gammaH*ha + p.chiH*ha*p.phi*(m1 + m2*p.phi2);
   Jac(8,3)=-p.chiH*ha*p.phi;
   Jac(8,4)=-p.chiH*ha*p.phi*p.phi2;
   Jac(8,5)=p.nu2*h;
   Jac(8,7)=p.nu2*c;
   Jac(8,8)=-p.gammaH - p.chiH*p.phi*(m1 + m2*p.phi2);
   Jac(9,5)=-p.rho(7)*p.r1*s;
   Jac(9,6)=p.r2*(1-s);
   Jac(9,9)=- c*p.rho(7)*p.r1 - p.r2*(g + 1);
   dy=dy';

   %Variational equation   
   dy(10:90)=Jac*Y;
   
end
