% This function outputs [dPhi/dA, dR/dA, dX/dA] 
function output = derivs1(a,phi,r,X,rd,Xd,Trh,Mphi,Gammaphi,GammaX,M,Md,xi,Bx,Bxd,sigmav,sigmavd)
global Mpl g Nbar gstar gd besselcutoff
denom=(phi + (r+rd)/a)^(1/2);
T=(30/(pi^2*gstar))^(1/4)*r^(1/4)/a*Trh;
Td=(30/(pi^2*gstar))^(1/4)*rd^(1/4)/a*Trh;
Ex=(M^2+9*T^2)^(1/2);
Xeq=(a/Trh)^3*g*T^3*(M/T)^2/(2*pi^2)*besselk(2,M/T);
if (M/T < besselcutoff)
    GammathX=GammaX/g*besselk(1,M/T)/besselk(2,M/T);
else
    GammathX=GammaX/g;
end
if (Md/Td< besselcutoff) 
GammathXinv= (M/Md)^2*GammaX/gd*besselk(1,M/Td)/besselk(2,Md/Td);
else
GammathXinv=(M/Md)^2*GammaX/gd*exp(-(M-Md)/Td)*(Md/M)^(1/2);    
end      
phidot=-(pi^2*gstar/30)^(1/2)*a^(1/2)*phi/denom;
% Note in lines below, we include Xeq (equil. abundance). 
% Uncomment lines below to take approx. Xeq = 0
%
rdot=(pi^2*gstar/30)^(1/2)*a^(3/2)*phi*(1 -(Md*Bxd + M*Bx)/Mphi)*(1 - xi)/denom + (3/(8*pi))^(1/2)*a^(-3/2)*sigmav*2*Ex*Mpl*(X^2 - Xeq^2)/denom;
%rdot=(pi^2*gstar/30)^(1/2)*a^(3/2)*phi*(1 -(Md*Bxd + M*Bx)/Mphi)*(1 - xi)/denom + (3/(8*pi))^(1/2)*a^(-3/2)*sigmav*2*Ex*Mpl*(X^2)/denom;
Xdot=-(3/(8*pi))^(1/2)*a^(-5/2)*sigmav*Mpl*Trh*(X^2 - Xeq^2)/denom + Gammaphi/Mphi*Bx*Nbar*(3/(8*pi))^(1/2)*(a)^(1/2)*Mpl/Trh*phi/denom + (3/(8*pi))^(1/2)*Mpl*a^(1/2)/Trh^2*(Xd*GammathXinv-X*GammathX)/denom;
%Xdot=-(3/(8*pi))^(1/2)*a^(-5/2)*sigmav*Mpl*Trh*(X^2)/denom + Gammaphi/Mphi*Bx*Nbar*(3/(8*pi))^(1/2)*(a)^(1/2)*Mpl/Trh*phi/denom + (3/(8*pi))^(1/2)*Mpl*a^(1/2)/Trh^2*(Xd*GammathXinv-X*GammathX)/denom;
output=[phidot;rdot;Xdot];
end