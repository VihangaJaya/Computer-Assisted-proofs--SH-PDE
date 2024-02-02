


x1 = solution_vector';
v = 1.1;

%solution_vector = [ 2.6923, -1.5635, 0.4999, 6.5938, -0.2978, 5.7463, 3.5358, 7.4444, -0.9423, 2.0971, -0.6555, 7.1535, -1.6423, 5.3296, 6.3223, 6.3323, -0.6843, 1.8764, 1.0559, 5.7208, 2.7051, 6.7076, 0.3931, 0.7970, 0.1790];

%Try 1.2
m =length(x1)-1;
%L = 100;
%alpha = 9.2735;

dx=ones(m+1,1);
iter=0;

while norm(dx,"inf")>1e-12 && iter<200
    f=function_handle(x1,m,L,alpha);
    df=jacobian(x1,L,alpha);
    dx=-df\f;
    x1=x1+dx;
    iter=iter+1;
end


mu = zeros(m+1,1);
for i = 0:m%changed the index
    mu(i+1)= -(((i)*pi)/L)^4 - 2*((i*pi)/L)^2 - (1-alpha);
end


Df = jacobian(x1,L,alpha);

%%
%z0
z0_vctor =zeros(m+1,1);
A_m = inv(jacobian(x1,L,alpha));
z0_b_vector = eye(m+1,m+1) - Df*A_m;


for i = 0:m
    sum_z0 = 0;
    for j = 0:m
        sum_z0 = sum_z0 + z0_b_vector(j+1,i+1);
    end
    z0_vctor(i+1) = sum_z0;
end

for k = 0:m
    z0_vctor(k+1) = (1/(v^k))*z0_vctor(k+1);
end
z0 = max(z0_vctor);

%%
%Z1 - rewritten


omega = zeros(3*m-1,1);
for n = 1: (3*m-1)
    if n > m+1
        omega(n) = v^(-n+1);
    end
end



convolution_vector = zeros(m+1,1);

for k = 0:m
    for k1 = -m:m
        for k2 = -m:m
            for k3 = m+2:3*m -1
                if (k1 +k2+k3) ==k
                    convolution_vector(k+1) = convolution_vector(k+1) + abs(x1(abs(k1)+1))*abs(x1(abs(k2)+1))*omega(k3);
                end
            end
        end
    end
end

product_matrix_conv = abs(A_m)*convolution_vector;

v_vector = v.*ones(m+1,1);
for k = 1:m+1
    v_vector(k) = v_vector(k)^(k-1);
end

z1 = 3*sum(abs(product_matrix_conv).*v_vector) + 48*v_norm(x1,v)^2/abs(mu(end));






%%

%z2

%calculating a_v wrong

z2_vector = zeros(m+1,1);
for i = 0:m
    sum_z2 = 0;
    for j = 0:m
        sum_z2 = sum_z2 + abs(A_m(j+1,i+1)) *v^L;
    end
    z2_vector(i+1) = sum_z2;
end

z2_vector = (1/v^m)*z2_vector;
k = max(z2_vector);
a_v = max(k,1/mu(m+1));

a_bar_norm = v_norm(A_m,v);
z2 = 96*a_v *a_bar_norm;

%%
%z3
z3 = 48*a_v;



%%

y_convolution1 = zeros(m+1,1);
y1 = 0;

F_m = function_handle(x1,m,L,alpha);
y1_convo = A_m * F_m;
for k = 0:m
    y1 = y1+ abs(y1_convo(k+1))*v^k;
end


y_convolution2 = zeros(2*m-2, 1);

for k = m:3*m-3
    U = -(((i)*pi)/L)^4 - 2*((i*pi)/L)^2 - (1-alpha);
    for k1 = -m:m
        for k2 = -m:m
            for k3 = -m:m
                if k1 + k2 + k3 == k
                    y_convolution2(k - m + 1) = y_convolution2(k - m + 1) + abs(x1(abs(k1) + 1)) * abs(x1(abs(k2) + 1)) * abs(x1(abs(k3) + 1));
                end
            end
        end
    end
    y2 = (1/abs(U)) * abs(y_convolution2(k - m + 1)) * v^k;
end


y = y1+y2;



%%


A=intval(A_m);
x=intval(x1);

Y0=y;

Z0=z0;

Z1=z1;

Z2=z2;
Z3 = z3;



Y0=intval(sup(Y0));
Z0=intval(sup(Z0));
Z1=intval(sup(Z1));
Z2=intval(sup(Z2));
Z3 = intval(sup(Z3));

rmin=NaN;
rmax=NaN;
if 1-Z0-Z1>0
    if 18*Z3*Z2*(1-Z0-Z1)*Y0-4*Z2^3*Y0+(Z2^2)*(1-Z0-Z1)^3-27*Z3^2*Y0^2 > 0
        intval_coefficients = [mid(Z3), mid(Z2), mid(Z0 +Z1 -1), mid(Y0)];
        roots_of_cubic = roots(intval_coefficients);
        two_largest = maxk(roots_of_cubic, 2);
        rmax_intval = intval(two_largest(1));
        rmin_interval = intval(two_largest(2));
        rmin = sup(rmin_interval);
        rmax = inf(rmax_intval);

        %rmax=inf(root);
        if rmin<rmax
            success=1;
            disp('success')
        else
            disp('failure: rmin > rmax')
        end
    else
        disp('failure: discriminant is negative')
    end
else
    disp('failure: linear term is positive')
end

r = [rmin, rmax];


%%

x = linspace(-10*pi,10*pi,1000);
plot(x,u(x,x1,L));



xlabel('Y'); % Replace 'X-axis Label' with your desired label for the x-axis
ylabel('U(Y)'); 


