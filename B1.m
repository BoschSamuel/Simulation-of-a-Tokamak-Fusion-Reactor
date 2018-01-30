function [Bx, By, Bz] = B1(input)
x = input(1);
y = input(2);
z = input(3);
distance = input(4);
a = input(5);
b = input(6);
I1 = input(7);
I2 = input(8);

mu = 4*pi*10^(-7);
r = sqrt(x^2+y^2); %distance from z-axis


%1st loop
z1 = z + distance/2;
k = sqrt(4*r*a/(z1^2+(a+r)^2));
[K,E] = ellipke(k);
Bz1 = mu*I1/(2*pi*sqrt(z1^2+(a+r)^2))*((a^2-z1^2-r^2)/(z1^2+(r-a)^2)*E+K);
Br1 = mu*z1*I1/(2*pi*r*sqrt(z1^2+(a+r)^2))*((z1^2+r^2+a^2)/(z1^2+(r-a)^2)*E-K);
Bx1 = Br1*x/r;
By1 = Br1*y/r;


%2nd loop
z2 = z - distance/2;
k = sqrt(4*r*a/(z2^2+(a+r)^2));
[K,E] = ellipke(k);
Bz2 = mu*I1/(2*pi*sqrt(z2^2+(a+r)^2))*((a^2-z2^2-r^2)/(z2^2+(r-a)^2)*E+K);
Br2 = mu*z2*I1/(2*pi*r*sqrt(z2^2+(a+r)^2))*((z2^2+r^2+a^2)/(z2^2+(r-a)^2)*E-K);
Bx2 = Br2*x/r;
By2 = Br2*y/r;

%central loop
z3 = z;
k = sqrt(4*r*b/(z3^2+(b+r)^2));
[K,E] = ellipke(k);
Bz3 = mu*I2/(2*pi*sqrt(z3^2+(b+r)^2))*((b^2-z3^2-r^2)/(z3^2+(r-b)^2)*E+K);
Br3 = mu*z3*I2/(2*pi*r*sqrt(z3^2+(b+r)^2))*((z3^2+r^2+b^2)/(z3^2+(r-b)^2)*E-K);
Bx3 = Br3*x/r;
By3 = Br3*y/r;


%total magnetic field - output
if (x==0)&&(y==0)
    Bx=0;
    By=0;
    Bz = Bz1 + Bz2 + Bz3;
else
    Bx = Bx1 + Bx2 + Bx3;
    By = By1 + By2 + By3;
    Bz = Bz1 + Bz2 + Bz3;
end



end