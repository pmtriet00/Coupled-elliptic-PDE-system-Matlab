
clear;


%Parameters according to page 16 - Dropbox Feb 2014 version
KA = 2.903;
KI = 0.018;
N = 6;
alpha = 2; %\alpha
Kr = 0.0033;
a0 = 0.33; 
H = 21;
tau = 0.8;
speed = 16.4;
z = 0.137;
r = 0.2;


R = 1280; %1/(tau*a0^H);
Kb = Kr*(1-a0)/a0;
p1 = alpha/(Kr+Kb);
q= Kr/(Kr+Kb);


V0top = N*H*(q-1)*R*q^H*(1+r*z+2*r*R*q^H)*speed^2;
V0bottom = (N*p1*q*(q-1) - (R*q^H+z)*(1+r*R*q^H))*(R*q^H+z)*(1 +r*R*q^H);
V0 = V0top / V0bottom;

K = V0*(KA - KI);



nx = 1500;
ny = 1500;
I1 = 2:nx/2+1;
I2 = nx/2+2:nx+1;
J = 2:ny;


xbeg = 0;
xend = 2;
x = linspace(xbeg, xend, nx+2);
x1 = linspace(xbeg, xbeg+ (xend - xbeg)/2, nx/2+1);

xi(I1) = 0.5*(x1(I1) + x1(I1-1));
xi(1) = x1(1) - 0.5*(x1(2) - x1(1));
xi(nx/2+1) = x1(nx/2+1) + 0.5*(x1(nx/2+1) - x1(nx/2));



S = nutrilinear(0,200,xbeg,xend/2,xi);
Smax = max(S(:));
ybeg = 0;
yend = 2; % ceil( (K*q/(1-q))^(1/N)*Smax) + 1;



tend =250;
dy= (yend - ybeg)/(ny+2);
dx= (xend - xbeg)/(nx+2);
dt = 10^-4*dy;





y = linspace(ybeg,yend,ny);



%Cell mesh for finite volume method
yi(J) = 0.5*(y(J) + y(J-1));
yi(1) = y(1) - 0.5*(y(2)-y(1));


%Density
p = zeros(nx/2+1,ny);
pnew = zeros(nx/2+1,ny);


u = zeros(nx+1,ny);
unew = zeros(nx+1,ny);


flux = zeros(nx+2,ny+1);


%Initialize u: delta y at y(0) = (K q /
%(1-q)^(1/N) * Smiddle

preadaptedy = 1; %(K*q / (1-q))^(1/N) * S(ceil(nx/4));


    for i = 1:nx
    for j = 1:ny
          
            if (0.999*preadaptedy <= y(j)) & (y(j) <= 1.001*preadaptedy)
         u(i,j) = 1;
        
            end
       
    end
    end
  

    
    
%Normalize mass to 1
u = u / (sum(u(:))*dx*dy);


%Recording the initial mean of u, use this for correction of mean later on
p(I1,J) = u(I1,J) + u(nx+2-I1,J);
px = zeros(1,nx/2+1);
for i = 1:nx/2+1
    
        px(i) = sum(p(i,:))*dy;   
    
end



% scaling z to 0:1
%    axis([0 xend/2 0 yend 0 1])

% scaling only x and y:
set(gca,'xlim',[0 xend/2])
set(gca,'ylim',[0 yend])

% for just xy view:
%    axis([0 xend/2 0 yend])

%Scaling z-axis
zmax = max(u(:));



f = zeros(nx+2,ny+1);

for i = 1 : nx/2+1
    for j = 2 : ny
   a(i,j) =  1/(1+K*((S(i) + KI)/((S(i) + KA) * y(j)))^N);
    end
end



for i = 1 : nx+1
    
        if i <= (nx / 2+1)
        f(i,J) = p1*y(J).*(q -  a(i,J)) ; %fplus
        else
        f(i,J) = p1*y(J).*(q -  a(nx+2-i,J)); %fminus
        end
      
    
end
    
lambda = zeros(nx/2+1,ny);

for i = 2 : nx/2+1
    lambda(i,J) = R*a(i,J).^H;
end


t = 0;
while t < tend
    
    unew(:,1) = 0;
    unew(:,ny) = 0;
    
      
    for i = 2 : nx+1
       for j = 2 : ny
            if f(i,j) > 0
        flux(i,j) = u(i,j-1)*f(i,j);
            else
        flux(i,j) = u(i,j)*f(i,j);
            end
       end
    end
    
     
       
    
     %Finite difference upwind scheme  
    unew(I1,J) = dt*[lambda(I1,J).*(u(nx+2-I1,J) - u(I1,J)) - speed*(u(I1,J) - u(I1-1,J))/dx - (flux(I1,J+1) - flux(I1,J))/dy ] + u(I1,J);        
    unew(I2,J) = dt*[lambda(nx+2-I2,J).*(u(nx+2-I2,J) - u(I2,J)) - speed*(u(I2,J) - u(I2-1,J))/dx - (flux(I2,J+1) - flux(I2,J))/dy ] + u(I2,J);   
    
    %Periodic boundary condition
    unew(1,:) = unew(nx+1,:);


%Verify mass is conserved
s = (sum(unew(:))- sum(unew(1,:)))*dx*dy; 


%Computing the density

 p(I1,J) = u(I1,J) + u(nx+2-I1,J);
 pnew(I1,J) = unew(I1,J) + unew(nx+2-I1,J);
 p(1,:) = u(1,:) + u(nx+1,:);
 pnew(1,:) = unew(1,:) + unew(nx+1,:);

 %Computing L-infty convergence
 dif = p - pnew;
pdif= max(abs((dif(:))));
 
t = t + dt;
fprintf('at time %10.2f out of %8.0f, the error is %12.10f, the mass is %1f \n',t,tend,pdif, s);
if (abs(t - 1) <= dt/2) | (abs(t - 50) <= dt/2)  | (abs(t - 100) <= dt/2) | (abs(t - 150) <= dt/2) | (abs(t - 200) <= dt/2) | (abs(t - 250) <= dt/2)
filename = [ 'myDataJuly22', num2str(t), '.mat' ];
save(filename);
 end
sfigure(1);
h=surf(x1,y,p','EdgeColor','none');       %plotting the velocity profile
    shading interp
    axis([xbeg xbeg+ (xend - xbeg)/2 0 yend 0 zmax*1.7])
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Chemical co-ordinate (y)')
    zlabel('Transport density profile (u) \rightarrow')
    drawnow;
    refreshdata(h)
    
 u = unew;

end

%Computing x-marginal density

px = zeros(1,nx/2+1);
for i = 1:nx/2+1
    
        px(i) = sum(p(i,:))*dy;   
    
end

sfigure(2);
plot(x1,px);
title('x-marginal density of bacteria');
