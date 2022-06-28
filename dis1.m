x0 = 0;
y0 = 2;
x1 = 1;
x2 = 1;
y1 = 3;
y2 = 1;

        
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
    reqAns = sqrt(x*x+y*y);
            
 elseif pr12_10 < 0 
    y = y0 - y1;
    x = x0 - x1;
    reqAns = sqrt(x*x+y*y);
            
 else
            
  x11 = x1x2;
  y11 = y1y2;
  x12 = x1x0;
  y12 = y1y0;
  mod = sqrt(x11*x11 + y11*y11);
  reqAns = abs(x11*y12 - y11*x12)/mod;
            
end
  
figure(1)        
plot(x0,y0,'rx')
hold on
plot(x1,y1,'rx')
hold on
plot(x2,y2,'rx')
axis equal