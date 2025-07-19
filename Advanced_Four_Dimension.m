%Hamed Baghal Ghaffari
%April 23 2023
%Construct Pseudo Zonal Spherical Monogenics


clear;
clc;

k=2;
syms x
Gegenoneovertwo=gegenbauerC(k,1,x);
diffGegenoneovertwo(x)=diff(Gegenoneovertwo,x);


Gegenthreeovertwo=gegenbauerC(k-1,2,x);
diffGegenthreeovertwo(x)=diff(Gegenthreeovertwo,x);

%a=zeros(3,k+1);
%for i=1:k+1
%x=rand(1,3);
%a(:,i)=x'/norm(x);
%end

%a=[-.962110407959442, -.211417294179026, -.172180982375314;-.0917140620969262, .989002178900472, -.11603112059399;.410386993978472, -.725295785362625, -.552746360308321;.57065831715715, .00126270293388618, .821186635675029;.195841697823354, -.0962555385036475, -.97590004647273;.140281073335604,-.592624038592789,.793169571668127;-.840938560236222,.131538146854321,.524900041713751];
%a=[-.962110407959442, -.211417294179026, -.172180982375314;-.0917140620969262, .989002178900472, -.11603112059399;.410386993978472, -.725295785362625, -.552746360308321;.57065831715715, .00126270293388618, .821186635675029;.195841697823354, -.0962555385036475, -.97590004647273];

% G=zeros(k+1,k+1);
% a=zeros((k+1)*(k+2)/2,4);
% 
% for i=1:(k+1)*(k+2)/2
%     theta1=(pi)*rand(1);
%     theta2=(pi)*rand(1);
%     phi=2*pi*rand(1);
%     x=sin(theta1)*cos(theta2);
%     y=sin(theta1)*sin(theta2)*cos(phi);
%     h=sin(theta1)*sin(theta2)*sin(phi);
%     z=cos(theta1);
%     a(i,1)=z;
%     a(i,2)=x;
%     a(i,3)=y;
%     a(i,4)=h;
% end

%prepared points for k=2
a=[-0.305395070204208  -0.367665390220912  -0.168257396130329  -0.862116825027788;
-0.391617999943390   0.037853646050709  -0.547734795122761   0.738369174473604; 
-0.169273440412214  -0.586050925623016  -0.791122360167905  -0.044902407394154;
0.605241395392630  -0.394272009905215  -0.420400404799795   0.549086455080301;
0.656068886383577  -0.304616197034322  -0.444683696465242  -0.528241421058805;
-0.938403153244574  -0.327094553087327  -0.110112132965520  -0.016852012787234];


%prepared points for k=3
%a=[0.268085164980490  -0.203278411976957  -0.248089276168423   0.908438188646723 ; -0.329785124910276   0.785481724563612  -0.511054346943372   0.114383942206812 ; -0.054461457024455   0.889809137772546   0.308762180393185   0.331571355810626 ; 0.827100090852648   0.100835312546690  -0.478931894353746   0.276336606380489 ; 0.267687635455224   0.345125905758713  -0.873446589164149  -0.215226612843613 ; 0.509774792458301   0.830859654652133  -0.219304281029268  -0.041322240593531 ; 0.321289766450786   0.648601452406885   0.514982994083370  -0.459218420487034 ; -0.788827525669338   0.335782566371205  -0.241475841876703  -0.454632401680188 ; -0.002447731204959   0.173623160389770   0.514142916926783   0.839944086092780 ; 0.931573491950828   0.178961197950981   0.242755801066834  -0.203010688809455];


% % for i=1:k+1
% %         for j=1:k+1
% %         G(i,j)=(2*k+1)*gegenbauerC(k,1/2,dot(a(i,:),a(j,:)))*;
% %         end
% % end


objectiveFun=0;

for l=1:(k+1)*(k+2)/2
for j=1:(k+1)*(k+2)/2
    if (l<j) 
    objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,2,dot(a(l,:),a(j,:))))^2;
    end
end
end


for p=1:2



for n=1:(k+1)*(k+2)/2



w=zeros(1,4);
for j=1:(k+1)*(k+2)/2
    if ~(j==n)  
    w=w-(((k+2)/2)^2*gegenbauerC(k,1,dot(a(n,:),a(j,:)))*diffGegenoneovertwo(dot(a(n,:),a(j,:))) ...
        -dot(a(n,:),a(j,:))*(gegenbauerC(k-1,2,dot(a(n,:),a(j,:))))^2 ...
        +(1-(dot(a(n,:),a(j,:)))^2)*gegenbauerC(k-1,2,dot(a(n,:),a(j,:)))*diffGegenthreeovertwo(dot(a(n,:),a(j,:)))).*(-a(j,:)+dot(a(n,:),a(j,:)).*a(n,:));
    end
end
%end


w=vpa(w/(norm(w)));

syms t
objectiveFun=0;
for j=1:(k+1)*(k+2)/2
    for l=1:(k+1)*(k+2)/2
        if ~(n==j) && ~(n==l) && (l<j)
            objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,2,dot(a(l,:),a(j,:))))^2;
        end
    end
end

for j=1:(k+1)*(k+2)/2
    if ~(n==j)
        objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2+(1-(dot(cos(t)*a(n,:)+sin(t)*w,a(j,:)))^2)*(gegenbauerC(k-1,2,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2;
    end
end

% 
% % Define the objective function as an anonymous function
f = @(t) double(subs(objectiveFun, t));
% 
% % Perform the optimization
[tmin, fval] = fminbnd(f, 0, 4*pi);

a(n,:)=cos(tmin)*a(n,:)+sin(tmin)*w;

objectiveFun=0;

for l=1:(k+1)*(k+2)/2
for j=1:(k+1)*(k+2)/2
    if (l<j) 
    objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,2,dot(a(l,:),a(j,:))))^2;
    end
end
end


end

n=randi([1,(k+1)*(k+2)/2]);
% 
w=zeros(1,4);
for j=1:(k+1)*(k+2)/2
    if ~(j==n)  
    w=w-(((k+2)/2)^2*gegenbauerC(k,1,dot(a(n,:),a(j,:)))*diffGegenoneovertwo(dot(a(n,:),a(j,:))) ...
        -dot(a(n,:),a(j,:))*(gegenbauerC(k-1,2,dot(a(n,:),a(j,:))))^2 ...
        +(1-(dot(a(n,:),a(j,:)))^2)*gegenbauerC(k-1,2,dot(a(n,:),a(j,:)))*diffGegenthreeovertwo(dot(a(n,:),a(j,:)))).*(-a(j,:)+dot(a(n,:),a(j,:)).*a(n,:));
    end
end
%end


w=vpa(w/(norm(w)));

syms t
objectiveFun=0;
for j=1:(k+1)*(k+2)/2
    for l=1:(k+1)*(k+2)/2
        if ~(n==j) && ~(n==l) && (l<j)
            objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,2,dot(a(l,:),a(j,:))))^2;
        end
    end
end

for j=1:(k+1)*(k+2)/2
    if ~(n==j)
        objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2+(1-(dot(cos(t)*a(n,:)+sin(t)*w,a(j,:)))^2)*(gegenbauerC(k-1,2,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2;
    end
end
% 
% % Define the objective function as an anonymous function
f = @(t) double(subs(objectiveFun, t));
% 
% % Perform the optimization
[tmin, fval] = fminbnd(f, 0, 4*pi);


a(n,:)=cos(tmin)*a(n,:)+sin(tmin)*w;

objectiveFun=0;

for l=1:(k+1)*(k+2)/2
for j=1:(k+1)*(k+2)/2
    if (l<j) 
    objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,2,dot(a(l,:),a(j,:))))^2;
    end
end
end


end
% 
% end
% 
% a
% 
a
objectiveFun=0;

for l=1:(k+1)*(k+2)/2
for j=1:(k+1)*(k+2)/2
    if (l<j) 
    objectiveFun=objectiveFun+((k+2)/2)^2*(gegenbauerC(k,1,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,2,dot(a(l,:),a(j,:))))^2;
    end
end
end

objectiveFun
syms e0 e12 e13 e14 e23 e24 e34

G=sym(zeros((k+1)*(k+2)/2));
for i=1:(k+1)*(k+2)/2
    G(i,i)=6.0*e0;
    for j=2:(k+1)*(k+2)/2
        if j>i
            G(i,j)=2*gegenbauerC(k,1,dot(a(i,:),a(j,:)))*e0+gegenbauerC(1,2,dot(a(i,:),a(j,:)))*((a(i,1)*a(j,2)-a(i,2)*a(j,1))*e12 + (a(i,1)*a(j,3)-a(i,3)*a(j,1))*e13+ ...
                    +(a(i,1)*a(j,4)-a(i,4)*a(j,1))*e14+(a(i,2)*a(j,3)-a(i,3)*a(j,2))*e23+(a(i,2)*a(j,4)-a(i,4)*a(j,2))*e24+(a(i,3)*a(j,4)-a(i,4)*a(j,3))*e34);
            G(j,i)=2*gegenbauerC(k,1,dot(a(i,:),a(j,:)))*e0+gegenbauerC(1,2,dot(a(i,:),a(j,:)))*(- (a(i,1)*a(j,2)-a(i,2)*a(j,1))*e12 - (a(i,1)*a(j,3)-a(i,3)*a(j,1))*e13 ...
                    - (a(i,1)*a(j,4)-a(i,4)*a(j,1))*e14 - (a(i,2)*a(j,3)-a(i,3)*a(j,2))*e23 - (a(i,2)*a(j,4)-a(i,4)*a(j,2))*e24 - (a(i,3)*a(j,4)-a(i,4)*a(j,3))*e34);
        end
    end
end
G=vpa(G,4)

syms e1 e2 e3 e12 e13 e23

T_mat_of_G=vpa(subs(G, {e12 e13 e14 e23 e24 e34}, {e1 e2 e3 e12 e13 e23}),4);

T_mat_of_G_plus=vpa(subs(T_mat_of_G, {e1 e2 e3}, {0 0 0}),4);

T_mat_of_G_minus=vpa(subs(T_mat_of_G,{e0 e12 e13 e23},{0 0 0 0}),4);

T_mat_of_G_tilde_plus=vpa(subs(T_mat_of_G_minus, {e1 e2 e3}, {e0 e12 e13}),4);

phi_T_mat_of_G_plus=vpa(subs(T_mat_of_G_plus, {e12 e13}, {-e12 -e13}),4);

phi_T_mat_of_G_tilde_plus=vpa(subs(T_mat_of_G_tilde_plus, {e12 e13}, {-e12 -e13}),4);

H=vpa([T_mat_of_G_plus,T_mat_of_G_tilde_plus;-phi_T_mat_of_G_tilde_plus,phi_T_mat_of_G_plus],4);

T_mat_of_H=vpa(subs(H, {e12 e13 e23}, {e1 e2 e12}),4);

T_mat_of_H_plus=vpa(subs(T_mat_of_H, {e1 e2}, {0 0}),4);

T_mat_of_H_minus=vpa(subs(T_mat_of_H,{e0 e12},{0 0}),4);

T_mat_of_H_tilde_plus=vpa(subs(T_mat_of_H_minus, {e1 e2}, {e0 e12}),4);

phi_T_mat_of_H_plus=vpa(subs(T_mat_of_H_plus, {'e12'}, {-e12}),4);

phi_T_mat_of_H_tilde_plus=vpa(subs(T_mat_of_H_tilde_plus, {'e12'}, {-e12}),4);

F=vpa([T_mat_of_H_plus,T_mat_of_H_tilde_plus;-phi_T_mat_of_H_tilde_plus,phi_T_mat_of_H_plus],4);

T_mat_of_F=vpa(subs(F, {e0 e12}, {1 1i}),4);

B=vpa((T_mat_of_F)^(-1/2),4);

T_mat_inverse_of_B=vpa(real(B)*e0+subs(B-real(B), {'1i'}, {e12}),4);


[r, c] = size(T_mat_inverse_of_B);

r2 = r / 2;
c2 = c / 2;

T_mat_inverse_of_B_plus = vpa(T_mat_inverse_of_B(1 : r2 , 1 : c2),4); 

T_mat_inverse_of_B_tilde_plus=vpa(T_mat_inverse_of_B(1 : r2 , c2 + 1 : c),4);

T_mat_inverse_of_B_tilde_plus_times_e1=vpa(expand(T_mat_inverse_of_B_tilde_plus*e1),4);

T_mat_inverse_of_B_minus=vpa(subs(T_mat_inverse_of_B_tilde_plus_times_e1, {e0 e1*e12}, {1 e2}),4);

Chi_inverse_T_mat_inverse_of_B=vpa(T_mat_inverse_of_B_minus+T_mat_inverse_of_B_plus,4);

T_mat_inverse_Chi_inverse_T_mat_inverse_of_B=vpa(subs(Chi_inverse_T_mat_inverse_of_B, {e1 e2 e12},{e12 e13 e23}),4);



r22 = r2 / 2;
c22 = c2 / 2;

T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_plus = vpa(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B(1 : r22 , 1 : c22)); 

T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_tilde_plus = vpa(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B(1 : r22 , c22 + 1 : c2),4);

T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_tilde_plus_timese1 = vpa(expand(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_tilde_plus*e1),4);

syms e123

T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_minus=vpa(subs(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_tilde_plus_timese1, {e0*e1 e1*e12 e1*e13 e1*e23}, {e1 e2 e3 e123}),4);


Chi_inverse_T_mat_inverse_Chi_inverse_T_mat_inverse_of_B=vpa(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_minus+T_mat_inverse_Chi_inverse_T_mat_inverse_of_B_plus,4);


syms e1234

T_mat_inver_Chi_inver_T_mat_inver_Chi_inver_T_mat_inver_B=vpa(subs(Chi_inverse_T_mat_inverse_Chi_inverse_T_mat_inverse_of_B, {'e1' 'e2' 'e3' 'e12' 'e13' 'e23' 'e123'}, {e12 e13 e14 e23 e24 e34 e1234}),4)





 