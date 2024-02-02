function vnorm = v_norm(a,v)
m = length(a)-1;
sum =0;
for i = 0:m
    sum =sum + abs(a(abs(i)+1)) * (v^abs(i));  
end
vnorm = sum ;
end
