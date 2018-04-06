
Finaloutputs=table2array(Finaloutputs);
Beta=Finaloutputs(:,1);
Ensemble= Finaloutputs(:,2);
pls=Finaloutputs(:,3);
rf=Finaloutputs(:,4);
svm=Finaloutputs(:,5);
vip_Mean=Finaloutputs(:,6);
vip_std=Finaloutputs(:,7);

upper=vip_Mean+vip_std;
lower=vip_Mean-vip_std;
bands=(1:1:500);

subplot(6,1,1)
stem(bands,Beta)
axis([0 500 0 1.1])
title('Ideal model','fontsize',20)
set(gca,'fontsize',14,'XTickLabel',[] )

subplot(6,1,2);
stem(bands, Ensemble)
axis([0 500 0 1.1])
title('Ensemble','fontsize',20)
set(gca,'fontsize',14,'XTickLabel',[])

subplot(6,1,3)
plot(bands,pls)
title('PLS','fontsize',20)
set(gca,'fontsize',14,'XTickLabel',[])

subplot(6,1,4)
plot(bands,rf)
title('RF','fontsize',20)
set(gca,'fontsize',14,'XTickLabel',[])

subplot(6,1,5)
plot(bands,svm)
title('SVM','fontsize',20)
set(gca,'fontsize',14,'XTickLabel',[])


subplot(6,1,6)
ax6=subplot(6,1,6)
plot(bands,vip_Mean)
ciplot(lower,upper,bands)
hold(ax6,'on')
plot(bands,ones(size(bands)),'LineWidth',2)
hold(ax6,'off')
title('PLSR_{MEAN}','fontsize',20)
xlabel('Predictors','FontSize', 30)
set(gca,'fontsize',14)
[ax,h2]=suplabel('Model coefficient','y'); 
set(h2,'FontSize',30) 














