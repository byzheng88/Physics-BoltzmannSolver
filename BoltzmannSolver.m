% This code uses the following notation:
% Phi = Modulus, R = visible sector radiation, X = visible sector LOSP
% Rd = dark sector radiation, Xd = dark sector dark matter
% See Readme for additional details 
%
% derivs1.m contain differential equations for Phi, R, X
% derivs2.m contain differential equations for Rd, Xd
% derivs1.m and derivs2.m must be included in same path as Boltzmannsolver
%
% --------------- Initialize Parameters-----------------------------------
% besselcutoff is value of x where K1(x), K2(x) are replaced with limiting
% expressions (see derivs1.m and derivs2.m)
global Mpl g Nbar gstar gd besselcutoff
Mpl=1.22*10^(19);logAmax=50;besselcutoff=500;imax=1;
%
% g, gd correspond to D.O.F. for X, Xd
g=2;gd=2;Nbar=1;
%
% In this code I have set gstar = gstard (gstar in the dark sector)
gstar=10.75;gstard=gstar;
%
% Bx is moduli branching ratio to LOSP, Bxd is branching ratio to LSP.
Bx=0.1;Bxd=0.1;
%
% xi (which should be between 0 and 1) parameterizes how much of the 
% modulus decays to dark radiation. 
% For xi = 0 no moduli decay to dark rad., for xi = 1 no moduli decay to
% vis rad.
xi=0.5;
% --------------- Initialize Masses (in GeV)-------------------------------
%Mphi = modulus mass
Mphi=50000;
%
% M is LOSP mass, Md is DM Mass
M=100; Md=0.1;
%
% --------------- Initialize Cross Sections and Decay Rates (in GeV)-------
% sigmav, sigmavd are thermal avg. ann. xsec; 
%sigmav for LOSP, sigmavd for LSP
sigmav=10^(-7);sigmavd=10^(-8);
% GammaX is X to Xd decay rate
GammaX=(10/137)^2*M/(8*pi);
% Trh (GeV) is temperature parameterizing modulus decay width
Trh =0.01;
Gammaphi=(4*pi^3*gstar/45)^(1/2)*Trh^2/Mpl;
% --------------- Set Intial Conditions-----------------------------------
% 
% Set Rmin to small value, or else will encounter issues in Bessel functions
Phimin=3/(8*pi)*Mpl^2*Mphi^2/Trh^4;Rmin=10^(-15);Xmin=0;
%
% ----------------Begin iterative ODE solving----------------------------- 
% For description of iterative ODE solving procedure, see Readme
%
callderiv1=@derivs1;
callderiv2=@derivs2;
% Set initial condition Rd = 0, Xd = 0
Rdfun=@(x) 0;
Xdfun=@(x) 0;
xdmax=0;
%
% jmax is max number of iterations in ODE solving 
jmax=4; 
% Iterate until answer converges within specified tolerance = tol
% If specified tolerance is unreached, stop after jmax iterations
tol=0.02;
% Solution is contained in arrays a, vis, loga_d, dark as follows:
% vis(i,1) is the value of Phi for logA = loga(i)
% vis(i,2) is the value of R for logA = loga(i)
% vis(i,3) is the value of X for logA = loga(i)
%
% dark(i,1) is the value of Rd for logA = loga_d(i)
% dark(i,2) is the value of Xd for logA = loga_d(i)
%
    for j=1:jmax,
    visderivs=@(loga,vis) callderiv1(exp(loga),vis(1),vis(2),vis(3),Rdfun(loga),Xdfun(loga),Trh,Mphi,Gammaphi,GammaX,M,Md,xi,Bx,Bxd,sigmav,sigmavd)*exp(loga);
%   ode15s is particularly good for stiff ODE's
    [loga,vis]=ode15s(visderivs,[0,logAmax],[Phimin Rmin Xmin]);
    phifun=@(x) interp1(loga,vis(:,1),x);
    rfun=@(x) interp1(loga,vis(:,2),x);
    Xfun=@(x) interp1(loga,vis(:,3),x);
    hidderivs=@(loga_d,dark) callderiv2(exp(loga_d),phifun(loga_d),rfun(loga_d),Xfun(loga_d),dark(1),dark(2),Trh,Mphi,Gammaphi,GammaX,M,Md,xi,Bx,Bxd,sigmav,sigmavd)*exp(loga_d);
    [loga_d,dark]=ode15s(hidderivs,[0,logAmax],[Rmin Xmin]);
    Rdfun=@(x) interp1(loga_d,dark(:,1),x);
    Xdfun=@(x) interp1(loga_d,dark(:,2),x);
        if and(j>=2,abs((dark(size(dark,1),2) - xdmax)/dark(size(dark,1),2))<tol)
            break
        end    
    xdmax=dark(size(dark,1),2);
    end
% End iterative ODE solving
%
% -------------- Plot X, Xd and Phi ------------------------------
% Plot as function of logA, normalized to max value
%
Xmax=max(vis(:,3));
Xdmax=max(dark(:,2));
Phimax=max(vis(:,1));
%
warning off
semilogy(loga,vis(:,3)/Xmax,'b',loga_d,dark(:,2)/Xdmax,'r',loga,vis(:,1)/Phimax,'k')
xlabel('logA')
ylabel('Fraction of Max. Value') 
legend('X','Xd','Phi')
axis([0 logAmax 10^(-10) 8])
warning on
%
%---------------------- Compute DM relic abundance ----------------------
% Obtain values at logAmax to rescale relic abundance to present day
%
Rf=vis(size(vis,1),2);
Tf=(30/(pi^2*gstar))^(1/4)*vis(size(vis,1),2).^(1/4)/exp(loga(size(loga,1)))*Trh;
Tdf=(30/(pi^2*gstar))^(1/4)*dark(size(dark,1),1).^(1/4)/exp(loga_d(size(loga_d,1)))*Trh;
Af=exp(loga(size(loga,1)));
Xdf= dark(size(dark,1),2);
Tnow=2.35*10^(-13);
%
% Dark matter relic abundance ($\Omega_{X_d} h^2$):
OmegaX=(((((Af/Trh*Md*Xdf/Rf)*Tf/Tnow))*4.17*10^(-5)));
S=sprintf('Fraction of observed dark matter: %s',OmegaX/0.12);
disp(S)