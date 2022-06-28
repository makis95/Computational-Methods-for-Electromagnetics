%random walk
clc;
clear all;
close all;

%open file and read the panels' coordinates
fileID = fopen('surfaces.txt','r');
[A , counts] = fscanf(fileID,'%f');
fclose(fileID);
fileID = fopen('gauss.txt','r');
[B , counts] = fscanf(fileID,'%f');
fclose(fileID);

% % surfaces
points = length(A)/3;
P0 = reshape(A',[3,points]);
P0 = P0';
%find number of 
change_surf = [0];
num_of_surf=1;


for i=1 : length(P0)
    if num_of_surf == P0(i,1)
        P(i,:) = P0(i,2:end);
        continue;
    else 
        num_of_surf = num_of_surf + 1; 
        change_surf = [change_surf i-1] ;
        P(i,:) = P0(i,2:end);
    end
end
change_surf = [change_surf i];


% %gauss surfaces
points = length(B)/3;
G0 = reshape(B',[3,points]);
G0 = G0';
%find number of 
change_gauss = [0];
num_of_gauss=1;

for i=1 : length(G0)
    if num_of_gauss == G0(i,1)
        G(i,:) = G0(i,2:end);
        continue;
    else 
        num_of_gauss = num_of_gauss + 1; 
        change_gauss = [change_gauss i-1] ;
        G(i,:) = G0(i,2:end);
    end
end
change_gauss = [change_gauss i];




% for i=1 : length(change_surf)-1
%    s = P(change_surf(i)+1:change_surf(i+1),:);
% 
%    figure(1);
%    line([s(:,1)],[s(:,2)])
%    hold on;
% end
%       
% for i=1 : length(change_gauss)-1
%    g = G(change_gauss(i)+1:change_gauss(i+1),:);
% 
% 
%    line([g(:,1)],[g(:,2)])
%    hold on;
% end
% axis equal;

%calculate the perimeter of each gauss surface
perimeter_gauss = zeros(1,(length(change_surf)-1));
k = 0 ;
for i = 1: length(change_gauss)-1
    % g is the gauss surface we are at
    g = G(change_gauss(i)+1:change_gauss(i+1),:);
    k = k + 1;
    for j = 1 :length(g)-1
       x1 = g(j,1);
       y1 = g(j,2);
       x2 = g(j+1,1);
       y2 = g(j+1,2);
       
       if x1 == x2 
           perimeter_gauss(k) = perimeter_gauss(k) + abs(y1-y2);
       elseif y1 == y2
           perimeter_gauss(k) = perimeter_gauss(k) + abs(x1-x2);
       end
       
    end
    
end

C = zeros(length(change_surf)-1);
Pij = zeros(length(change_surf)-1);
Srfc = 0;
E = 8.85*10^(-12);
Nij = zeros(length(change_surf)-1);
for i = 1: length(change_gauss)-1
    % g is the gauss surface we are at as a matrix
    gaussian = G(change_gauss(i)+1:change_gauss(i+1),:)

    Srfc = Srfc + 1 %number of surface
    
    stop = 0;
    N=0; % walk counter of every surface 
    while stop == 0
        
        Pi = E*perimeter_gauss(Srfc);
       
        %a random point of the surface
        [p,nxy] = random_point(gaussian);
        [SQR,position,d] = square_create(p,P,change_surf);
        L = SQR(2,1) - SQR(1,1);
        perimeterSQR = abs(4*L); %perimeter of the 1st square
        Pi = Pi * perimeterSQR;

        
        %a random point at the first square
        [p1,nxy] = random_point(SQR); %random point from the square
        [SQR1,position1,d]  = square_create(p1,P,change_surf); %maximum square of x point
        L1 = SQR1(2,1) - SQR1(1,1);
        area1 = L1*L1; %area of the square
        
        %%chech if the point belogns to a surface
        while (area1 < 0.5) | area1 > 300
            [p1,nxy] = random_point(SQR); %random point from the square
            [SQR1,position,d] = square_create(p1,P,change_surf); %maximum square of x point
            L1 = SQR1(2,1) - SQR1(1,1);
            area1 = L1*L1; %area of the square
        end
         
        perimeterSQR1 = abs(4*L1);
         
        gax = 0;
        gay = 0;
        for n = 1:1:100
             gax = gax + ((2*n*pi*cos((n*pi)/2)*sinh((n*pi)/2)*sin(n*pi*p1(1)/L))/(L^2*sinh(n*pi)));
             gay = gay + ((2*n*pi*sin((n*pi)/2)*cosh((n*pi)/2)*sin(n*pi*p1(1)/L))/(L^2*sinh(n*pi)));
        end
        Pi = Pi * ( -nxy(1)*gax - nxy(2)*gay);
        
         %%% ** in the algorithm
        end_of_walk = 1;
        while end_of_walk ~= 0
           

            Pi = Pi*perimeterSQR1;
            
            [p2,nxy] = random_point(SQR1); %random point from the square
            [SQR2,position,d]  = square_create(p2,P,change_surf); %maximum square of x1 point
            L2 = SQR2(2,1) - SQR2(1,1);
            area2 = L2*L2; %area of the square
            perimeterSQR2 = abs(4*L2);%perimeter of square
            
            while  area2 > 300
                [p2,nxy] = random_point(SQR1); %random point from the square
                [SQR2,position,d] = square_create(p1,P,change_surf); %maximum square of x point
                L2 = SQR2(2,1) - SQR2(1,1);
                area2 = L2*L2; %area of the square
                perimeterSQR2 = abs(4*L2);%perimeter of square
            end
            
           gax = 0;
           for n = 1:1:100
                gax = gax + ((2*cos((n*pi)/2)*sinh((n*pi)/2)*sin(n*pi*p2(1)/L1))/(L1*sinh(n*pi)));
           end
           Pi = Pi * gax;
           
           % % check if we are at a surface
           if d >= 10^(-15) %we check the distanse of the random point from a surface
               end_of_walk = 1;
               
               % for the plot
               SQR0 = SQR1;
               
               %for the next iteration
               SQR1 = SQR2;
               perimeterSQR1 =  perimeterSQR2;
               L1 = L2;

           else %we hit surface

               % % chech at which surface we hit
               k = 0;
               for counter = 1 : length(change_surf)-1
                    k = k + 1;
                    if  change_gauss(counter)+1 <= position & position <= change_gauss(counter+1)
                        j0 = k;
                    end
               end
               
               Nij(Srfc,j0) = Nij(Srfc,j0) + 1;
               N = N + 1;
               Pij(Srfc,j0) = Pij(Srfc,j0) + Pi;
               C_prev = C(Srfc,j0);
               C(Srfc,j0) = Pij(Srfc,j0)/Nij(Srfc,j0);
                  
               if N>2000 %stop the walks from a surface after 10000 walks
                    stop = 1;
               else
                    stop = 0;
               end
         
               end_of_walk = 0;
               
           end
           
% %            % plots 
% %            figure(Srfc)
% %            for z=1 : length(change_surf)-1
% %                s = P(change_surf(z)+1:change_surf(z+1),:);
% %                line([s(:,1)],[s(:,2)])
% %                hold on;
% %            end
% % 
% %            for z=1 : length(change_gauss)-1
% %                g = G(change_gauss(z)+1:change_gauss(z+1),:);
% %                line([g(:,1)],[g(:,2)],'Color','cyan','LineStyle','--')
% %                hold on;
% %            end
% % 
% %             plot(p(1),p(2),'ro')
% %             hold on;
% %             plot(p1(1),p1(2),'bx')
% %             hold on;
% %             plot(p2(1),p2(2),'bx')
% %             hold on;
% %             plot(SQR(:,1),SQR(:,2));
% %             hold on;
% %             plot(SQR0(:,1),SQR0(:,2));
% %             hold on;
% %             plot(SQR1(:,1),SQR1(:,2));
% %             hold on;
% %             plot(SQR2(:,1),SQR2(:,2));
% %             hold on;
% %             axis equal
% %             
        end % 2nd while

    end %1st while
end



% A function that take the matrix of the surface 
% Nx2 and returns the p=(p(1),p(2)) = (x,y) which is the 
% the random point of the surface
function [p nxy] = random_point(g)

    g0 = randi([1 length(g)-1]);
    
    % vectors nx and ny
    if g0 == 1 
        nxy = [-1 0];
    elseif g0 == 2
        nxy = [0 1];
    elseif g0 == 3
        nxy = [1 0];
    else
        nxy = [0 -1];
    end
       
        
    x0 = g(g0,1);
    y0 = g(g0,2);
    x1 = g(g0+1,1);
    y1 = g(g0+1,2);
    if x0==x1

        if y1>y0
            p_y = y0+(y1-y0).*rand(1);
            p = [x0 p_y];
        else
            p_y = y1+(y0-y1).*rand(1);
            p = [x0 p_y];
        end
        
    else

        if x1>x0
            p_x = x0+(x1-x0).*rand(1);
            p = [p_x y0];
        else
            p_x = x1+(x0-x1).*rand(1);
            p = [p_x y0];
        end
    end

end


% A function that creates the squares
function [square pos d] = square_create(p,surface,change_surf)

x0 = p(1);
y0 = p(2);
d1 = 100000;
position = 0;
n = 0;
k = 0;
for i=1 : length(change_surf)-1
%     position = position + 1;
    s0 = surface(change_surf(i)+1:change_surf(i+1),:);
    
    for j = 1:length(s0)
        n = n+1;
        if j == length(s0);
            x1 = s0(j,1);
            x2 = s0(1,1);
            y1 = s0(j,2);
            y2 = s0(1,2);
        else
            x1 = s0(j,1);
            x2 = s0(j+1,1);
            y1 = s0(j,2);
            y2 = s0(j+1,2);
        end
        
        x1x2 = x2 - x1;
        y1y2 = y2 - y1;
        
        x2x0 = x0 - x2;
        y2y0 = y0 - y2;
        
        x1x0 = x0 - x1;
        y1y0 = y0 - y1;
        
        pr12_20 = x1x2*x2x0 + y1y2*y2y0;
        pr12_10 = x1x2*x1x0 + y1y2*y1y0;
        
        reqAns = 0;
        
        if pr12_20 > 0 
            y = y0 - y2;
            x = x0 - x2;
            d0(n) = 1000;
            pos(n) = n;
        elseif pr12_10 < 0 
            y = y0 - y1;
            x = x0 - x1;
           d0(n) = 1000;
           pos(n) = n;

        else
 
            x11 = x1x2;
            y11 = y1y2;
            x12 = x1x0;
            y12 = y1y0;
            mod = sqrt(x11*x11 + y11*y11);
            d0(n) = abs(x11*y12 - y11*x12)/mod;
            pos(n) = n;

        end
        
    end

end
 
    
% find at which corner is the random point near
dmin=10000;
for i = 1: length(surface)
    d1(i) = max(abs(x0- surface(i,1)),abs(y0-surface(i,2)));
    pos1(i) = i;
end

[min1,i1] = min(d0);
[min2,i2] = min(d1);

if min1 > min2
    d = min2;
    pos = i2;
else
    d = min1;
    pos = i1;
end


sqr1(1,:) = [(p(1) - d) (p(2) - d)];
sqr1(2,:) = [(p(1) + d) (p(2) - d)];
sqr1(3,:) = [(p(1) + d) (p(2) + d)];
sqr1(4,:) = [(p(1) - d) (p(2) + d)];
sqr1(5,:) = [(p(1) - d) (p(2) - d)];

square = sqr1;
end


