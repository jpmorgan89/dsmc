function freemolecule
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all
N = 20000;%no of molecules
M = 20;%no of cells

beta_u=1;%beta_U/beta_L
%beta_l=1
u_u=0;%upper wall u-velocity 

edges = linspace(0,1,M+1);

delta_t = (1/M)/2/(1/sqrt(pi));

z=rand([N 1]);

[u, v, w]=maxwell(N);%maxwell-distributed vels

%p=1;
%p=0.5;
p=0; %probability of diffuse reflection

endtime=4/(1/sqrt(pi));%average velocity is 1/sqrt(pi)

bottomhit=[];
ztrack=[];
for i=1:endtime/delta_t
    for j = 1:N
        z(j) = z(j)+delta_t*w(j);
        
        if z(j)>1
            time_rem=(z(j)-1)/w(j);
            if rand<p
                %diffuse reflection
                [unew,vnew,~]=maxwell(1);
                wnew = -(beta_u^-1)*sqrt(-log(1-rand));
                u(j)=beta_u^-1.*unew+u_u;
                v(j)=beta_u-1.*vnew;
                w(j)=wnew;
            else
                %specular reflection
                w(j)=-w(j);
            end
            z(j)=1+w(j)*time_rem;
        end
        
        if z(j)<0
            time_rem=(z(j))/w(j);
            bottomhit = [bottomhit w(j)];%keeping track of bottom wall hits
            
            if rand<p
                [unew,vnew,~]=maxwell(1);
                wnew = sqrt(-log(1-rand));
                u(j)=unew;
                v(j)=vnew;
                w(j)=wnew;
            else
                w(j)=-w(j);
            end
            z(j)=w(j)*time_rem;
        end
        
        
    end
    ztrack = [ztrack z(1)];
end

subplot(4,1,1)
binsize=.1;
hist(w(z<edges(2)),-5.05:binsize:5.05)
hold on
title('Cell Next to Top Wall')
ylabel('Frequency')
xlabel('w/(2RT_L)^{1/2}')
plot(linspace(-5.05,5.05,1000),sum(z<edges(2))*binsize*pi^(-1/2)*exp(-linspace(-5.05,5.05,1000).^2),'g')

subplot(4,1,2)
binsize=.1;
hist(w(z>edges(end-1)),-5.05:binsize:5.05)
hold on
title('Cell Next to Bottom Wall')
ylabel('Frequency')
xlabel('w/(2RT_L)^{1/2}')
plot(linspace(-5.05,5.05,1000),sum(z>edges(end-1))*binsize*pi^(-1/2)*exp(-linspace(-5.05,5.05,1000).^2),'g')

subplot(4,1,3)
binsize=.1;
hist(w(z>edges(floor(M/2))&z<edges(1+floor(M/2))),-5.05:binsize:5.05)
hold on
title('Cell In Middle')
ylabel('Frequency')
xlabel('w/(2RT_L)^{1/2}')
plot(linspace(-5.05,5.05,1000),sum(z>edges(floor(M/2))&z<edges(1+floor(M/2)))*binsize*pi^(-1/2)*exp(-linspace(-5.05,5.05,1000).^2),'g')

subplot(4,1,4)
binsize=.1;
hist(bottomhit,-5.05:binsize:5.05)
hold on
title('Bottom wall strikes')
ylabel('Frequency')
xlabel('w/(2RT_L)^{1/2}')
plot(linspace(-5.05,5.05,1000),-length(bottomhit)*binsize*2*linspace(-5.05,5.05,1000).*exp(-linspace(-5.05,5.05,1000).^2),'g')

figure

hist(z,.025:.05:.975)

figure

plot(ztrack)
end
% useless comment