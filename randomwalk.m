clc
close all
clear
fileID = fopen('conductor1.txt','r');
[A,count1] = fscanf(fileID, '%d' )
fclose(fileID);
fileID = fopen('gauss1.txt','r');
[A1,count12] = fscanf(fileID, '%d' )
fclose(fileID);
fileID = fopen('conductor2.txt','r');
[B,count2] = fscanf(fileID, '%d' )
fclose(fileID);
fileID = fopen('gauss2.txt','r');
[B1,count22] = fscanf(fileID, '%d' )
fclose(fileID);
fileID = fopen('conductor3.txt','r');
[C,count3] = fscanf(fileID, '%d' )
fclose(fileID);
fileID = fopen('gauss3.txt','r');
[C1,count33] = fscanf(fileID, '%d' )
fclose(fileID);
ep=exp(-20);
A=A';
A1=A1';
x1=[];
y1=[];
x1_g=[];
y1_g=[];
x2_g=[];
y2_g=[];
x3_g=[];
y3_g=[];
for i=1:2:count1-1 
    x1 =[x1 A(i+1)];
    y1=[y1 A(i)];
end
for i=1:2:count12-1 
    x1_g =[x1_g A1(i+1)];
    y1_g=[y1_g A1(i)];
end

for i=1:2:count22-1 
    x2_g =[x2_g B1(i+1)];
    y2_g=[y2_g B1(i)];
end
for i=1:2:count33-1 
    x3_g =[x3_g C1(i+1)];
    y3_g=[y3_g C1(i)];
end
plot (x1,y1);
hold on
plot (x1_g,y1_g);
hold on
plot (x2_g,y2_g);
hold on
plot (x3_g,y3_g);

B=B';
x2=[];
y2=[];

for i=1:2:count2-1 
    x2 =[x2 B(i+1)];
    y2=[y2 B(i)];
end
hold on
plot (x2,y2);
C=C';
x3=[];
y3=[];

for i=1:2:count3-1 
    x3 =[x3 C(i+1)];
    y3=[y3 C(i)];
end


hold on
plot (x3,y3);

l1= length (x1_g);
l2= length (x2_g);
l3= length (x3_g);


per1 = perimeter(x1_g,y1_g);
per2=perimeter(x2_g,y2_g);
per3=perimeter(x3_g,y3_g);
E=8.85*(10^-12);
k=1;
P22=0;
P21=0;
P23=0;
N2=0;
c22=[];
c21=[];
c23=[];
l=1
m=1;
i=1;
tol=0.001;
while k<50
[p,p1,ra]=create_random_gauss(l2,x2_g,y2_g,tol);
[P,con1,con2,con3]=random_wlk(per2,x1,x2,x3,y1,y2,y3,E,ep,p,p1,tol)
if con1
    l=l+1;
    P21=P21+P;
    N2=N2+1;
        c21(l)=P21/N2

end
if con2
    m=m+1;
    P22=P22+P;
    N2=N2+1;
    c22(m)=P22/N2;

end
if con3
    i=i+1;
    P23=P23+P;
    N2=N2+1;
    c23(i)=P23/N2;
end
k=k+1;
end

function [P,con1,con2,con3]=random_wlk(per,x1,x2,x3,y1,y2,y3,E,ep,p,p1,tol)
l1c= length (x1);
l2c= length (x2);
l3c= length (x3);
P=per*E;
min_distance=min_dist(x1,x2,x3,y1,y2,y3,p,p1);
[r,s]=create_rect(min_distance,p);
P=P*s;
[p,p1,ra]=create_random_gauss(5,r(:,1),r(:,2),tol);
min_distance=min_dist(x1,x2,x3,y1,y2,y3,p,p1);

 while (min_distance<0.01)
      [p,p1,ra]=create_random_gauss(5,r(:,1),r(:,2),tol);
      min_distance=min_dist(x1,x2,x3,y1,y2,y3,p,p1);
 end


 gax=0;
 a=abs(r(1,1)-r(2,1))
[nx,ny,x]= check(ra,p)

 for n=1:1:101
     gax=gax+ ((2*n*pi*cos((n*pi)/2)*sinh((n*pi)/2)*sin((n*pi*x)/a))/(a*a*sinh(n*pi)));
 end
 n=0;
 gay=0;
 for n=1:1:101
     gay=gay+ ((2*n*pi*sin((n*pi)/2)*cosh((n*pi)/2)*sin((n*pi*x)/a))/(a*a*sinh(n*pi)));
 end
 P=P*((-nx*gax)-(ny*gay));
 n=0;
ga=0;
k=1;
i=1;
check1=0;
  check2=0;
  check3=0;
   con1=0;
   con2=0;
   con3=0;
[r,s]=create_rect(min_distance,p);
P=P*s;
 while min_distance>ep  ;
 [p,p1,ra]=create_random_gauss(5,r(:,1),r(:,2),tol);
  min_distance=min_dist(x1,x2,x3,y1,y2,y3,p,p1);
  if(min_distance<0.01 && min_distance>ep) || (min_distance>10 );
  while (min_distance<0.01 && min_distance>ep) || (min_distance>10 ); 
      [p,p1,ra]=create_random_gauss(5,r(:,1),r(:,2),tol);
      min_distance=min_dist(x1,x2,x3,y1,y2,y3,p,p1);
  end
  end
  
   [nx,ny,x]=check(ra,p);
   a=abs(r(1,1)-r(2,1))
  for n=1:1:101
     ga=ga+ ((2*sin((n*pi)/2)*sinh((n*pi)/2)*sin((n*pi*x)/a))/(a*sinh(n*pi)));
  end
  P=P*ga;
  if min_distance>ep
   [r,s]=create_rect(min_distance,p); 
   P=P*s;
  end
   k=k+1;
   
 end

 if min_distance<ep
  p(1)=round (p(1),3);
  p(2)=round (p(2),3);
    for i=1:l1c-1
        [check1(i),en]=checkPointOnSegment([x1(i) y1(i)],[x1(i+1) y1(i+1)],p);
        if check1(i)
            con1=1;
        end
    end
    for i=1:l2c-1
        [check2(i),en]=checkPointOnSegment([x2(i) y2(i)],[x2(i+1) y2(i+1)],p);
        if check2(i)
            con2=1;
        end
    end
    for i=1:l3c-1
        [check3(i),en]=checkPointOnSegment([x3(i) y3(i)],[x3(i+1) y3(i+1)],p);
    if check3(i)
            con3=1;
        end
    end
    
  end


i=1;
end

function p=perimeter(x,y)
l=length(x)
p=0;
for i=1:l-1
    p=p+abs(x(i+1)-x(i))+abs(y(i+1)-y(i))
    i
end 
end

function [checkPt, onEnd] = checkPointOnSegment(A,B,C)



checkPt = false;
onEnd = false;


AB = B - A;
AC = C - A;


if cross(AB, AC) == 0
    % calculate the dotproduct of (AB, AC) and (AB, AB) to see point is now
    % on the segment
    dotAB = dot(AB, AB);
    dotAC = dot(AB, AC);
    % on end points of segment
    if dotAC == 0 || dotAC == dotAB
        onEnd = true;
        checkPt = true;
    % on segment
    elseif dotAC > 0 && dotAC < dotAB
        checkPt = true;
    end
end
end
function z = cross(a, b)
    z = a(1)*b(2) - a(2)*b(1);
end


 
 function [squart,s]=create_rect(min_distance,p) 
       width = min_distance; % whatever
height = min_distance; % whatever...
xCenter = p(1); % Wherever...
yCenter = p(2); % Wherever...
xLeft = xCenter - width;
yBottom = yCenter - height;
r=rectangle('Position', [xLeft, yBottom, width*2, height*2], 'EdgeColor', 'b');
squart(1,:)=[xLeft yBottom];
squart(2,:)=[xLeft+width*2 yBottom];
squart(3,:)=[xLeft yBottom+width*2];
squart(4,:)=[xLeft+width*2 yBottom+width*2];
squart(5,:)=[xLeft yBottom];
s= 2*abs(squart(1,1)-squart(2,1))+2*abs(squart(1,2)-squart(3,2))
grid on;
 end
 
 function [nx,ny,x]=check(ra,p)
  if ra==1
     nx=-1;
     ny=0;
     x=p(1);
 end
  if ra==2
     nx=0;
     ny=1;
      x=p(2);
     
  end
      if ra==3
     nx=1;
     ny=0;
      x=p(1);
      end
      if ra==4
     nx=0;
     ny=-1;
     x=p(2);
      end
     if ra==5
        nx=1;
        ny=0;
         x=p(1);
     end
 end
 
 function min_distance=min_dist(x1,x2,x3,y1,y2,y3,p,p1)
      l1= length (x1);
      l2= length (x2);
      l3= length (x3);
   
      z1=0;
      z2=0;
      z3=0;
      for i=1:l1-1
       z1=[0 z1];
      end
      for i=1:l2-1
       z2=[0 z2];
      end
      for i=1:l3-1
       z3=[0 z3];
      end
      
           q1 =[x1', y1'
                x2',y2'
                x3',y3'
                ];
            q2 =[x1', y1',z1'
                x2',y2',z2'
                x3',y3',z3'
                ];
          
          for i=1:l1-1
           [distance(i),vert(i)] = point_to_line_seg(p,q1(i,:), q1(i+1,:));
           distance1(i) = point_to_line(p1,q2(i,:), q2(i+1,:));
          end 
          for i=l1+1:l2+l1-1
           [distance(i-1),vert(i-1)] = point_to_line_seg(p,q1(i,:), q1(i+1,:));
           distance1(i-1) = point_to_line(p1,q2(i,:), q2(i+1,:));
          end 
          for i=l2+l1+1:l3+l1+l2-1
           [distance(i-2),vert(i-2)] = point_to_line_seg(p,q1(i,:), q1(i+1,:));
           distance1(i-2) = point_to_line(p1,q2(i,:), q2(i+1,:));
          end 
          
           j=1;
           
          flag=0;
       for i=1:l1+l2+l3-3
           if vert(i)>=0
              dist(j)=distance(i);
              j=j+1;
              flag=1;
           end
       end
       if flag
            min_distance=min(dist);
       end
       
      
       
       
       for i=1:l1-1
           if distance1(i)>distance1(i+1)
               dst(i)=distance1(i);
           else
               dst (i)=distance1(i+1);
           end
           if distance1(l1-1)>distance1(1)
               dst(l1-1)=distance1(l1-1);
           else 
               dst(l1-1)=distance1(1);
           end
       end
       for i=l1+1:l1+l2-1
           if distance1(i-1)>distance1(i)
               dst(i-1)=distance1(i-1);
           else
               dst (i-1)=distance1(i);
           end
           if distance1(l1+l2-2)>distance1(l1)
               dst(l1+l2-2)=distance1(l1+l2-2);
           else 
               dst(l1+l2-2)=distance1(l1);
           end
       end
           for i=l1+l2:l1+l2+l3-3
           if distance1(i-1)>distance1(i)
               dst(i-1)=distance1(i-1);
           else
               dst (i-1)=distance1(i);
           end
           if distance1(l1+l2+l3-3)>distance1(l1+l2-1)
               dst(l1+l2+l3-3)=distance1(l1+l2+l3-3);
           else 
               dst(l1+l2+l3-3)=distance1(l1+l2-1);
           end
           end
            min_dst=min(dst);
            
            if flag
            min_distance=min (min_dst,min_distance);
            else 
                min_distance=min_dst;
                flag=0;
            end
 end
function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end
function [dist,vert] = point_to_line_seg(x, a, b)
d_ab = norm(a-b);
d_ax = norm(a-x);
d_bx = norm(b-x);
vert=dot(a-b,x-b)*dot(b-a,x-a);
if dot(a-b,x-b)*dot(b-a,x-a)>=0
    A = [a,1;b,1;x,1];
    dist = abs(det(A))/d_ab;        
else
    dist = min(d_ax, d_bx);
end         
end

function [p,p1,ra]=create_random_gauss(l,x1_g,y1_g,tol)

ra=randi (l);

random1=[x1_g(ra) y1_g(ra)];
if ra==l
random1n=[x1_g(2) y1_g(2)];
else 
    random1n=[x1_g(ra+1) y1_g(ra+1)];
end

if random1(1)==random1n(1)
    if random1(2)<random1n(2)
        p(1)=random1(1);
        a = random1(2);
b = random1n(2);
p(2)= (b-a)*rand(1) + a;
        
    else
         p(1)=random1(1)
   b = random1(2);
a = random1n(2);
p(2)= (b-a)*rand(1) + a;
    end
else
    if random1(1)<random1n(1)
        a = random1(1);
b = random1n(1);
p(1)= (b-a)*rand(1) + a;
        p(2)=random1(2);
    else 
   b = random1(1);
a = random1n(1);
p(1)= (b-a)*rand(1) + a;
    p(2)=random1(2);
    end
end
%for i=1:l-1
 %  if p(1)==x1_g(i) && p(2)==y1_g(i)
  %      if x1_g(i)<x1_g(i+1)
   %     p(1)==p(1)+tol;
    %    end
    %end
%end
plot (p(1),p(2),'*')
p1=[p(1) p(2) 0];
end