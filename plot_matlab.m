close all
clear all

re = load('pearsonr','-ascii'); 
se = load('slope','-ascii');
fisher = load('fisher','-ascii');


for i = 1:size(re,1)
    f(i) = re(i,1);
    ravg(i) = re(i,2);
    rmax(i) = re(i,3);
    rmin(i) = re(i,4);
    rvaravg(i) = re(i,5);
    rvarmax(i) = re(i,6);
    rvarmin(i) = re(i,7);
end

for i = 1:size(se,1)
    fs(i) = se(i,1);
    savg(i) = se(i,2);
    smax(i) = se(i,3);
    smin(i) = se(i,4);
    svaravg(i) = se(i,5);
    svarmax(i) = se(i,6);
    svarmin(i) = se(i,7);
end

for i = 1:size(fisher,1)
    ff(i) = fisher(i,1);
    zavg(i) = fisher(i,2);
    zvar(i) = fisher(i,3);
    zarea(i) = fisher(i,4);
end


errorbar(100*(f),ravg,(rmax-ravg),(ravg-rmin),'ko-','Markersize',8) % 'Markersize',20)
xlabel('Sample percent of population (S_p = 100n/N)','FontSize',14,'FontWeight','bold')
ylabel('Difference in Pearson R expectation','FontSize',14,'FontWeight','bold')
pause
errorbar(100*(fs),ravg,(smax-savg),(savg-smin),'ko-','Markersize',8) % 'Markersize',20)
ylabel('Difference in slope expectation','FontSize',14,'FontWeight','bold')
pause
errorbar(100*(f),100*rvaravg,100*(rvarmax-rvaravg),100*(rvaravg-rvarmin),'ko-','Markersize',8) % 'Markersize',20)
ylabel('Percent relative difference in Pearson R variance','FontSize',12,'FontWeight','bold')
pause
errorbar(100*(fs),100*svaravg,100*(svarmax-svaravg),100*(svaravg-svarmin),'ko-','Markersize',8) % 'Markersize',20)
ylabel('Percent relative difference in slope variance','FontSize',12,'FontWeight','bold')
pause
plot(100*ff,zavg,'ko-','Markersize',8)
ylabel('Difference in Fisher transformed expectation','FontSize',12,'FontWeight','bold')
pause
plot(100*ff,100*zvar,'ko-','Markersize',8)
ylabel('Percent relative difference in Fisher variance','FontSize',12,'FontWeight','bold')
pause
plot(100*ff,zarea,'ko-','Markersize',8)
ylabel('Error in normal approximation','FontSize',12,'FontWeight','bold')











