function u = u(x,x1,L)

u = 0;
for i = 0:size(x1,1)-1

    u = u + x1(i+1)*cos(i*x*pi/L);

end
