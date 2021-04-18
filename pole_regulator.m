ku = 1/(4.187);
kv = 1;
VR = 2;
Q = 1;
L = VR / Q;
V1 = 10;
V2 = 12;
s=tf('s') ;  % specificerar att överföringsfunktionen är i Laplace
G=ku*exp(-L*s)*1/(((V1*V2)/(Q*Q))*s*s + (((V1 + V2)/Q)*s) + 1);  % överföringsfunktionen G(s)

s = tf('s');
GP = G;
h = 1;
HP=c2d(GP,h);
fprintf('\nNegative notation for HP')
HP.variable = 'z^-1'
fprintf('\nPositive notation for HP')
HP.variable = 'z'
[B,A]=tfdata(HP,'v')
pol = 0.6;
P =poly([pol pol]);
Kr = polyval(P,1) / polyval(B,1)

% lös ekvationssystemet
syms c1 c2 c3 d0 d1
eqn1 = c1 + A(2) == P(2);
eqn2 = c2 + A(2)*c1 + A(3) == P(3);
eqn3 = c3 + A(2)*c2 + A(3)*c1 + B(2)*d0 == 0;
eqn4 = A(3)*c2 + A(2)*c3 + B(3)*d0 + B(2)*d1 == 0;
eqn5 = A(3)*c3 +B(3)*d1 == 0;
[M,Yp] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5], [c1, c2, c3, d0, d1]);
X = linsolve(M,Yp);

% lägg till värderna som decimaltal
c_1 = double(X(1))
c_2 = double(X(2))
c_3 = double(X(3))
d_0 = double(X(4))
d_1 = double(X(5))

Htot.variable = 'z'

%Process' step response.
figure('name','Step Response of Process','numbertitle','off')
step(GP)

D = [d_0  d_1 ];
D_block = tf( D,[1 0],h);
fprintf('\nPositive notation for D block')
D_block.variable = 'z'

C = [1 c_1 c_2 c_3 ];
C_block = tf([1 0 0 0], C, h);
fprintf('\nPositive notation for C block')
C_block.variable = 'z'

Kr = polyval(P,1) / polyval(B,1)
[GpNum,GpDen]=tfdata(GP,'v')

N = 100;
r1 = [zeros(1,4),ones(1,N-4)];
y1 = zeros(1, 100);
u1 = zeros(1, 100);
y2 = zeros(1, 100);
u2 = zeros(1, 100);
y3 = zeros(1, 100);

for k = 5:N
    
    u1(k) = -P(2)*u1(k-1) - P(3)*u1(k-2) + Kr + A(2)*Kr + A(3)*Kr;
    y1(k) = B(2)*u1(k-3) + B(3)*u1(k-4) - A(2)*y1(k-1) - A(3)*y1(k-2);

end
k1 = 10:1:(N+9);
plot(k1,y1,k1,u1)
legend('y(k)','u(k)')
title('2 Poler i Z = 0,6 och resten i origo')
xlabel('k') 
ylabel('y(k) / u(k)') 