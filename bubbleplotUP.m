function h=bubbleplotUP(x,y,z,color,sf,color_flag,N,maxxim,Til)
%function h=bubbleplot(x,y,z,color,sf)
%PURPOSE
%  Scatter plot in x and y using arbitrary-color "." as
%  plot symbol - whose point-by-point size is proportional
%  the the magnitude of z
%INPUT
%  x     - n-dimension vector
%  y     - n-dimension vector
%  z     - n-dimension vector (used to size the plot symbols)
%  color - plot symbol color (must be a 3-element vector
%          with elements in range 0 1) ................... default: alternating 10 colors
%  sf    - plot symbol size scale factor ................. default: 35 

n = nargin;
if n<6 | isempty(color_flag), color_flag = 0; end
if n<5 | isempty(sf), sf = 35; end
if n<4 | isempty(color),
   myco=[ 0     0     1.00
          0     0.50  0
          1.00  0     0
          0     0.75  0.75
          0.75  0.75  0.75
          0.75  0     0.75
          0.25  1.00  0.25
          0.75  0.75  0
          0.25  0.25  0.25
          0.50  0.50  0.50 ];
else   
   col = color(:)';
   myco = [col;col;col;col;col;col;col;col;col;col];
end
% scf=sf./(max(z(:))-min(z(:)));
% zp=round( (z-min(z(:)) ) .*scf +8  ); 
zp=3*z;
I = find(isnan(zp));
zp(I) = eps; 
if color_flag==1
ZP=Til-maxxim+1;
% zp=4*ones(1,length(zp));
xx=linspace(0,1,N);
yy=linspace(1,0,N);
zz=linspace(1,1,N);
dd=[xx' yy' zz'];
Agg=dd(ZP,:);
% zp=4*ones(1,length(zp));
end

for i = 1:length(x)
   cc = i;
   if cc>10
      cc = mod(cc,10)+1;
   end
   if (y(i)~=1) 
   plot(x(i),y(i),'.','MarkerSize',zp(i),'Color',myco(cc,:));
   hold on;
   elseif y(i)==1 
      if color_flag==1 
       % plot(x(i),y(i),'.','MarkerSize',zp(i),'Color',Agg);
       plot(x(i),y(i),'.','MarkerSize',11.55,'Color',Agg);
       hold on;
      else
       plot(x(i),y(i),'.','MarkerSize',zp(i),'Color',myco(cc,:));
       hold on;
      end
   end
end
hold off
axis('tight')
ax=axis ;             % L R B T
xr=.03*(ax(2)-ax(1));
yr=.03*(ax(4)-ax(3));% scale axis based on data range
axis([ax(1)-xr ax(2)+2*xr ax(3)-yr ax(4)+2*yr]);
hold off;
h = gca;