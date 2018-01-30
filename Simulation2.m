clear

initial_position = [2.3 0 0]; % x z y
initial_speed = [-0.1 -0.15 +0]; % vx vy vz

q=1.60217*10^(-19); %charge
m=1.67262*10^(-27); % mass
c = 3*10^8;
I_coils = 80; %current in the coil
N = 15000; %number of windings
I_coils = N*I_coils;
I_plasma = 1.0e6;
inside_tokamak = true;
g = 9.81; 

a = 1.5; %radius of each coil
b = 0.8; %radius of central region

t_final = 1.0e-6; %duration of sim.
dt = 2.0e-11; %step size

disp(['Speed = ', num2str((initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)*100) , '% speed of light'])
disp(['Energy = ', num2str((1/sqrt(1-initial_speed(1)^2-initial_speed(2)^2-initial_speed(3)^2)-1)*m*c^2/q/10^6), ' MeV'])

x = initial_position(1);
y = initial_position(2);
z = initial_position(3);
vx = initial_speed(1);
vy = initial_speed(2);
vz = initial_speed(3);

t = (0:dt:t_final);

plot3(x,y,z,'.b');
xlim([-4 4])
ylim([-4 4])
zlim([-3.5 3.5])
xlabel('x');
ylabel('y');
zlabel('z');
hold on
grid on

for i = 0:15
    plotCircle3D([(a+b)*cos(pi/16+i*pi/8) (a+b)*sin(pi/16+i*pi/8) 0], [-sin(pi/16+i*pi/8) cos(pi/16+i*pi/8) 0], a, true);
end


gamma = 1/sqrt(1-(vx^2+vy^2+vz^2));
vx = vx*c;
vy = vy*c;
vz = vz*c;
count = 0;

waitforbuttonpress

for i = 0:t_final/dt
    
    phi =  atan2(y,x);
    distance = sqrt( z^2 + (x-(a+b)*cos(phi))^2 + (y-(a+b)*sin(phi))^2 );
    if (distance>a)||(x^2+y^2>(b+2*a)^2)||(x^2+y^2<b^2)
        disp('The particle has escaped!');
        break
    end
    
    [Bx, By, Bz] = B2([x y z a b I_coils I_plasma]);  %magnetic field strength calc.
    %Bx=0; Bz=0; By=0.0000002;   %uniform magnetic field
    
    ax = q/(gamma*m)*(vy*Bz-vz*By);
    ay = q/(gamma*m)*(vz*Bx-vx*Bz);
    az = q/(gamma*m)*(vx*By-vy*Bx)-g;
    
    vx = vx + ax*dt;
    vy = vy + ay*dt;
    vz = vz + az*dt;
    
    %sim. error correction
    vx = vx*sqrt(initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)/sqrt((vx^2+vy^2+vz^2)/c^2);
    vy = vy*sqrt(initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)/sqrt((vx^2+vy^2+vz^2)/c^2);
    vz = vz*sqrt(initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)/sqrt((vx^2+vy^2+vz^2)/c^2);
    
    x = x + vx*dt;
    y = y + vy*dt;
    z = z + vz*dt;
    count = count+1;
    
    if mod(count,60)==0
        plot3(x,y,z,'.b', 'linewidth',1000);
    end
    
    if mod(count,180)==0
        drawnow
    end

end


hold off