%% Set default parameter values according to Table 1 of the manuscript

function p=setBaselineParams()

p.nu=0.1;
p.nu2=0.01;
p.phi=0.01;
p.phi2=10;
p.chiH=1;

% Platelet scaling parameter
p.rho=1; 

% Decay rates
p.gammaA=1;
p.gammaM=0.01;
p.gammaM2=0.01;
p.gammaG=1;
p.gammaH=0.1;
p.gammaE=0.1;

% Rate parameters
p.kn=0.01;
p.kg=0.1;
p.kh=0.1;
p.km=0.0001;
p.km1=30;
p.km2=0.3;
p.r1=1;
p.r2=1;

% Saturation constants
p.betaA=0.1;
p.betaN=0.1;
p.betaC=0.12;
p.betaG=0.01;

end