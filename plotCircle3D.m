function plotCircle3D(center,normal,radius, grayline)
%This function needs three input parameters: 
%*Circle radius 
%*position of the circle midpoint 
%*vector perpendicular to the surface in which the circle is plotted. 
%With help of scatter3, the circle is plotted into the current axes and the handle for the new plot is returned to the user


theta=0:0.02:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'r-', 'linewidth',6);
if grayline
    plot3(points(1,:),points(2,:),points(3,:).*0-3.5, 'color',[0.5 0.5 0.5], 'linewidth',1);
end

end