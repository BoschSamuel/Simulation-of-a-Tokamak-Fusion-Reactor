function [Bx, By, Bz] = B2(input)
x = input(1);
y = input(2);
z = input(3);
a = 2*input(4);
b = input(5);
I_coils = input(6);
I_plasma = input(7);

mu = 4*pi*10^(-7);

Bx=0; By=0; Bz=0;

%magnetic field of the coils
for i = 0:15
    theta = pi/16+i*pi/8; %angle between the i-th coil and the x-axis
    
    if abs(sin(theta))>0.01
        r1_ = x/(cos(theta)-sin(theta)*tan(atan2(y,x)-theta))-a-b; 
        r1 = sqrt(r1_^2+z^2);
        z1 = ((b+a)*cos(theta)+r1_*cos(theta)-x)/sin(theta);
    else
        r1_ = x-b-a;
        r1 = sqrt(r1_^2+z^2);
        z1 = y;
    end
    
    k = sqrt(4*r1*a/(z1^2+(a+r1)^2));

    [K,E] = ellipke(k);

    Bz1_ = mu*I_coils/(2*pi*sqrt(z1^2+(a+r1)^2))*((a^2-z1^2-r1^2)/(z1^2+(r1-a)^2)*E+K); %Bz1_ is the magnetic field in the coil frame
    Br1_ = mu*z1*I_coils/(2*pi*r1*sqrt(z1^2+(a+r1)^2))*((z1^2+r1^2+a^2)/(z1^2+(r1-a)^2)*E-K); %Br1_ is the magnetic field in the cail frame

    Bx1 = -sin(theta)*Bz1_+Br1_*r1_/r1*cos(theta); %normal coordinates
    By1 = cos(theta)*Bz1_ + sin(theta)*Br1_*r1_/r1; %normal coordinates
    Bz1 = Br1_*z/r1; %normal coordinates

    %adding the field of a single coil to the total field
    Bx = Bx+Bx1;
    By = By+By1;
    Bz = Bz+Bz1;
    
    if abs(Bx)<(5e-12) %JIK.
        Bx=0;
    end
    if abs(By)<(5e-12) %JIK.
        By=0;
    end
    if abs(Bz)<(5e-12) %JIK.
        Bz=0;
    end
    
end


%magnetic field of the plasma current
sigma = a/3; %parameter of the Gauss curve
phi =  atan2(y,x);
distance = sqrt( z^2 + (x-(a+b)*cos(phi))^2 + (y-(a+b)*sin(phi))^2 ); %distance to centre of plasma ring
I2_r_plasma = I_plasma*erf(distance/(sigma*sqrt(2)));

r = sqrt(x^2+y^2);
k = sqrt(4*r*(a+b)/(z^2+((a+b)+r)^2));
[K,E] = ellipke(k);
Bz_plasma = mu*I2_r_plasma/(2*pi*sqrt(z^2+((a+b)+r)^2))*(((a+b)^2-z^2-r^2)/(z^2+(r-(a+b))^2)*E+K);
Br_plasma = mu*z*I2_r_plasma/(2*pi*r*sqrt(z^2+(b+r)^2))*((z^2+r^2+(a+b)^2)/(z^2+(r-(a+b))^2)*E-K);
Bx_plasma = Br_plasma*x/r;
By_plasma = Br_plasma*y/r;

if (distance>0.0001) %JIK.
    Bx = Bx+Bx_plasma;
    By = By+By_plasma;
    Bz = Bz+Bz_plasma;
end



end
