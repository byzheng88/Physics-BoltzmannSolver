% This function outputs [dRd/dA, dXd/dA] 
function output = derivs2(a,phi,r,X,rd,Xd,Trh,Mphi,Gammaphi,GammaX,M,Md,xi,Bx,Bxd,sigmav,sigmavd)
global Mpl g Nbar gstar gd besselcutoff
denom=(phi + (r+rd)/a)^(1/2);
T=(30/(pi^2*gstar))^(1/4)*r^(1/4)/a*Trh;
Td=(30/(pi^2*gstar))^(1/4)*rd^(1/4)/a*Trh;
Ex=(M^2+9*T^2)^(1/2);
Exd=(Md^2+9*Td^2)^(1/2);
Xdeq=(a/Trh)^3*gd*Td^3*(Md/Td)^2/(2*pi^2)*besselk(2,Md/Td);
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
% Note in lines below, we include Xdeq (equil. abundance). 
% Uncomment lines below to take approx. Xeq = 0
%
rddot=(pi^2*gstar/30)^(1/2)*a^(3/2)*phi*(1 -(Md*Bxd + M*Bx)/Mphi)*xi/denom + (3/(8*pi))^(1/2)*a^(-3/2)*sigmavd*2*Exd*Mpl*(Xd^2 - Xdeq^2)/denom+(3/(8*pi))^(1/2)*a^(3/2)*Mpl/(Trh^3*denom)*(Ex - Exd)*(GammathX*X - GammathXinv*Xd);
%rddot=(pi^2*gstar/30)^(1/2)*a^(3/2)*phi*(1 -(Md*Bxd + M*Bx)/Mphi)*xi/denom + (3/(8*pi))^(1/2)*a^(-3/2)*sigmavd*2*Exd*Mpl*(Xd^2)/denom+(3/(8*pi))^(1/2)*a^(3/2)*Mpl/(Trh^3*denom)*(Ex - Exd)*(GammathX*X - GammathXinv*Xd);
Xddot=-(3/(8*pi))^(1/2)*a^(-5/2)*sigmavd*Mpl*Trh*(Xd^2 - Xdeq^2)/denom + Gammaphi/Mphi*Bxd*Nbar*(3/(8*pi))^(1/2)*(a)^(1/2)*Mpl/Trh*phi/denom -(3/(8*pi))^(1/2)*Mpl*a^(1/2)/Trh^2*(Xd*GammathXinv-X*GammathX)/denom;
%Xddot=-(3/(8*pi))^(1/2)*a^(-5/2)*sigmavd*Mpl*Trh*(Xd^2)/denom + Gammaphi/Mphi*Bxd*Nbar*(3/(8*pi))^(1/2)*(a)^(1/2)*Mpl/Trh*phi/denom -(3/(8*pi))^(1/2)*Mpl*a^(1/2)/Trh^2*(Xd*GammathXinv-X*GammathX)/denom;
output=[rddot;Xddot];