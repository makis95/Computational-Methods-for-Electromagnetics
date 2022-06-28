%method of moments

clc;
clear all;
close all;

%open file and read the panels' coordinates
fileID = fopen('panels.txt','r');
[A , counts] = fscanf(fileID,'%f');
fclose(fileID);



%define the number of panels
panels_count = counts/9; 
voltage = 1;

%P is the matrix of the coordinates of panels
%where the first number is the number of the panel we are
P0 = reshape(A',[9,panels_count]);
P0 = P0';
change_panel = [0];
pnl=1;
P = [];
for i=1 : length(P0)
    if pnl == P0(i,1)
        P(i,:) = P0(i,2:end);
        continue;
    else 
        pnl = pnl + 1; 
        change_panel = [change_panel i-1] ;
        P(i,:) = P0(i,2:end);
    end
end
change_panel = [change_panel i];


x = [];
y = [];
for i=1:panels_count
    for j=1:2:8
        x = [x P(i,j)];
        y = [y P(i,j+1)];
    end
end

N = 2; %number of subpanels we want to create

%create subpanels
counter = 1;
for i=1 : 4 : length(x) 
   
   x_step = (x(i+2) - x(i))/N ;
   y_step = (y(i+1) - y(i))/N ;
   x_start = x(i);
   y_start = y(i);
   
   for k = 1:N^2
        
       xi(counter) = x_start;
       xi(counter+1) = x_start; 
       xi(counter+2) = x_start + x_step;
       xi(counter+3) = x_start + x_step;
       
       
       yi(counter) = y_start;
       yi(counter+1) = y_start + y_step; 
       yi(counter+2) = y_start;
       yi(counter+3) = y_start + y_step;
       y_start = y_start + y_step;
       
       counter = counter+4;
       
       if mod(k,N)==0
           x_start = x_start + x_step;
           y_start = y(i);
       end
   end
   
%    for j = x(i): x_step : x(i+2)
%  
%        for k = y(i) : y_step : y(i+1)
% 
%             xi(counter) = j;
%             yi(counter) = k;
%             counter = counter +1;
%             
%        end
%    end
end


%find centers
j=1;
for i=1:4:length(xi)-2
    center(j,1) = (xi(i+2) + xi(i))/2;
    center(j,2) = (yi(i+1)+ yi(i))/2;
    j = j+1;
end

figure(2);
plot(xi,yi,'xr','Linewidth',1);
hold on;
plot(center(:,1),center(:,2),'ob','Linewidth',1);


%surface of panels and hx - hy for every square
dS = [];
j = 1;
for i=1:4:length(xi)
    dS(j) = (xi(i+3)-xi(i))*(yi(i+3)-yi(i));
    hxi(j) = xi(i+3)-xi(i);
    hyi(j) = yi(i+3)-yi(i);
    j = j+1;
end
change_panel = N^2*change_panel;

%for loop for each surface
C = [];
for i=1:length(change_panel)-1
   
   %voltage initialize depending on the surface we are 
   vi = zeros(1,length(dS));
   for n=change_panel(i)+1:change_panel(i+1)
      vi(n) = 1;
   end
   
   
  %for loop for each panel  
  for j=1:length(dS)
      
      %create equation of Ö(rj)
      equ = [];

      for n=1:length(dS)
        if(j ~= n)
            equ(n) = 1/(4*pi*(8.85*(10^(-12))))*(dS(n)/(norm([center(j,1) center(j,2)] - [center(n,1) center(n,2)])));
           
        else
            hx = hxi(j);
            hy = hyi(j);
            equ(n) = 2/(4*pi*(8.85*(10^(-12))))*(hx*log((hy+sqrt(hx^2+hy^2))/hx)+hy*log((hx+sqrt(hx^2+hy^2))/hy));            
        end
        
      end %equation loop
      % %create matrix the every panel equation
      M(j,:) = equ;
      
  end %panel loop
  
  % %solve the system for the surface
  surf = M\vi';
  
  % % find Cij
  for j=1: length(change_panel)-1
      from = floor(change_panel(j));
      to = floor(change_panel(j+1));
      C(i,j) = surf(from+1:to)'*dS(from+1:to)';
  end
    
  
end %surface loop

