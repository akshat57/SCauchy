
del=1.0006e-3;
ut=4.9968e-2;
tnu=del/ut; 

cau=load('omegamn.dat'); 
dsplus=(cau(:,1)-cau(1,1))/tnu; 
caup=cau*tnu;
plot(dsplus,caup(:,2),'-r','LineWidth',2)
hold on
plot(dsplus,caup(:,3),'-b','LineWidth',2)
plot(dsplus,caup(:,4),'-g','LineWidth',2)
xlabel('$\delta s^+$','FontSize',20,'interpreter','latex') 
ylabel('${\bf E}\left[\tilde{\omega}_i^+(s)\right]$','FontSize',20,'interpreter','latex') 
hl=legend('$i=x$','$i=y$','$i=z$'); 
set(hl,'interpreter','latex','FontSize',15,'Location','SouthEast') 
print('omegamn.jpg','-djpeg') 


var=load('varomega.dat');
varp=var*tnu^2; 
figure; semilogy(dsplus,varp(:,1),'-r','LineWidth',2)
hold on
semilogy(dsplus,varp(:,2),'-b','LineWidth',2)
semilogy(dsplus,varp(:,3),'-g','LineWidth',2)
xlabel('$\delta s^+$','FontSize',20,'interpreter','latex') 
ylabel('${\rm Var}\left[\tilde{\omega}_i^+(s)\right]$','FontSize',20,'interpreter','latex') 
hl=legend('$i=x$','$i=y$','$i=z$'); 
set(hl,'interpreter','latex','FontSize',15,'Location','SouthEast') 
print('varomega.jpg','-djpeg') 
