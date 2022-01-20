close all
clear all
clc
tic
n=input('Enter step of qw: ')
a=[[ 0.     0.     0.     0.     0.     0.     0.     0.707  0.     0....
   0.     0.     0.     0.     0.    -0.707]
 [ 0.     0.     0.     0.     0.     0.     0.707  0.     0.     0....
   0.     0.     0.     0.    -0.707  0.   ]
 [ 0.     0.     0.     0.     0.707  0.     0.     0.     0.     0....
   0.     0.    -0.707  0.     0.     0.   ]
 [ 0.     0.     0.     0.     0.     0.707  0.     0.     0.     0....
   0.     0.     0.    -0.707  0.     0.   ]
 [ 0.707  0.     0.     0.     0.     0.     0.     0.    -0.707  0....
   0.     0.     0.     0.     0.     0.   ]
 [ 0.     0.707  0.     0.     0.     0.     0.     0.     0.    -0.707...
   0.     0.     0.     0.     0.     0.   ]
 [ 0.     0.     0.707  0.     0.     0.     0.     0.     0.     0....
  -0.707  0.     0.     0.     0.     0.   ]
 [ 0.     0.     0.     0.707  0.     0.     0.     0.     0.     0....
   0.    -0.707  0.     0.     0.     0.   ]
 [ 0.     0.     0.     0.     0.707  0.     0.     0.     0.     0....
   0.     0.     0.707  0.     0.     0.   ]
 [ 0.     0.     0.     0.     0.     0.707  0.     0.     0.     0....
   0.     0.     0.     0.707  0.     0.   ]
 [ 0.     0.     0.     0.     0.     0.     0.707  0.     0.     0....
   0.     0.     0.     0.     0.707  0.   ]
 [ 0.     0.     0.     0.     0.     0.     0.     0.707  0.     0....
   0.     0.     0.     0.     0.     0.707]
 [ 0.     0.     0.707  0.     0.     0.     0.     0.     0.     0....
   0.707  0.     0.     0.     0.     0.   ]
 [ 0.     0.     0.     0.707  0.     0.     0.     0.     0.     0....
   0.     0.707  0.     0.     0.     0.   ]
 [ 0.     0.707  0.     0.     0.     0.     0.     0.     0.     0.707...
   0.     0.     0.     0.     0.     0.   ]
 [ 0.707  0.     0.     0.     0.     0.     0.     0.     0.707  0....
   0.     0.     0.     0.     0.     0.   ]];
qa=[1 0];
qb=[1 0];
qc=[1 0];
syms x y
h=[cos((x)/2) (exp(1*i*(y))*sin((x)/2))];
k=kron(h,qa);
kk=kron(k,qb);
j=kron(kk,qc);
syms jj(x,y)
jj(x,y)=transpose(j)
x=0:0.1:6.28;
y=0:0.1:6.28;
qq=length(x)
qc=length(y)
z=zeros(qq,qc);
for tx=1:qq
    for ty=1:qc
        ll=single((jj(x(tx),y(ty))));
        ll=(a^n)*ll;
        gg=(ll).*conj(ll);
        sd=std(gg);
        z(tx,ty)=sd;   
    end
end
minval = min(z(:))
[p q]=find(z==minval)
min_theta=x(p)
min_phi=y(q)
%meshc(x,y,z)
surfc(x,y,z)
%hidden off
colormap('hsv')
xlabel('x')
ylabel('y')
zlabel('z')
toc
if (min_theta>3.14) 
    min_theta=min_theta-3.14
end
best_ugate=[cos((min_theta)/2)  -1*sin((min_theta)/2);(exp(1*i*(min_phi)/2))*(sin((min_theta)/2)) (exp(1*i*(min_phi)/2))*(cos((min_theta)/2))]
initialized_qcoin=best_ugate*[1;0]




%axis([0 4 0 8 0.08 0.126])
%subplot(2,1,1)
%meshc(x,y,z)
%subplot(2,1,2)
%fsurf(w, [0 3.14,0 6.28])
