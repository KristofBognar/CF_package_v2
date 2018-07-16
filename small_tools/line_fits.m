function [intercept,slop,slop_nlm,mdl_lm,mdl_nlm] = line_fits(x,y)

mdl_lm = fitlm(x,y);

myfun = @(b,x)b(1)*x;
beta =[1];
mdl_nlm = fitnlm(x,y,myfun,beta);
figure;hold all;
plot(x,y,'.');
plot(x,predict(mdl_lm,x),'color',[0.9 0.1 0.5]);
plot(x,predict(mdl_nlm,x),'color',[0.1 0.1 0.9]);
plot(x,x,'k--');
intercept = mdl_lm.Coefficients.Estimate(1);
slop = mdl_lm.Coefficients.Estimate(2);
slop_nlm = mdl_nlm.Coefficients.Estimate(1);
textbp(['y = ' num2str(slop) '*x + ' num2str(intercept)],'color',[0.9 0.1 0.5]);
textbp(['y = ' num2str(slop_nlm) '*x'],'color',[0.1 0.1 0.9]);
textbp(['N = ' num2str(numel(x)) ]);
textbp(['R = ' num2str(corr(x,y)) ]);