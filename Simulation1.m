clear

initial_position = [0 0 -.8]; % x z y
initial_speed = [0.8 -0.48 0.3595]; % vx vy vz


q=1.60217*10^(-19); %charge
m=9.10938215*10^(-31); % mass
c = 3*10^8; %speed of light
I1 = 80; %current in the solenoid
N1 = 1500; %number of windings
I1 = N1*I1;

I2 = 100; %current in the central solenoid
N2 = 8000; %number of windings
I2 = N2*I2;

distance = 10; %distance between solenoids
a = 1; %radius of each coil
b = 4; %radius of central coil

t_final = 3.0*10^(-7); %duration of sim.
dt=1.5*10^(-12); %step size

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
xlim([-5 5])
ylim([-5 5])
zlim([-5.2 5.2])
xlabel('x');
ylabel('y');
zlabel('z');
hold on
grid on

plotCircle3D([0 0 -distance/2], [0 0 1], a, false);
plotCircle3D([0 0 distance/2], [0 0 1], a, false);
plotCircle3D([0 0 0], [0 0 1], b, false);

gamma = 1/sqrt(1-(vx^2+vy^2+vz^2));
vx = vx*c;
vy = vy*c;
vz = vz*c;
count = 0;

waitforbuttonpress

for i = 0:t_final/dt
    
    [Bx, By, Bz] = B1([x y z distance a b I1 I2]);  %magnetic field strength calc.
    %Bx=0; Bz=0; By=0.0000002;   %uniform magnetic field
    
    ax = q/(gamma*m)*(vy*Bz-vz*By);
    ay = q/(gamma*m)*(vz*Bx-vx*Bz);
    az = q/(gamma*m)*(vx*By-vy*Bx);
    
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
    
    if mod(count,200)==0
        plot3(x,y,z,'.b', 'linewidth',1000);
    end
    
    if mod(count,400)==0
        drawnow
    end

end

hold off
