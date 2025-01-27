FS='FontSize';
x1 = -1.65:.02:1.65; y1 = -1.5:.02:1.5;
[xx,yy] = meshgrid(x1,y1); ss = xx+1i*yy;
u = @(s) -1 + 0.5*real((s+1).*log(s+1)-(s-1).*log(s-1));
hold off
cval = linspace(-1,0.5,11);
contour(x1,y1,u(ss),cval,'k','linewidth',1.4)
hold on
contour(x1,y1,u(ss),(-1+log(2))*[1 1],'r','linewidth',1.4)
%set(gca,'xtick',-2:.5:2,'ytick',-.5:.5:.5), grid on
ylim([-1.5 1.5]), axis equal, grid on
hold on, plot(.5255i,'ok')
text(0.05,.63,'0.52552491457i')
title(['Runge region for equispaced interpolation ' ...
       'in the limit n -> \infty'],FS,9)
hold off
print -dpdf 'equipot_uniform.pdf'

figure()
u = @(s) log(abs(s+1i*sign(imag(s)).*sqrt(1-s.^2))) - log(2);
xgrid = -1.5:.02:1.5; ygrid = -0.91:.02:0.91;
[xx,yy] = meshgrid(xgrid,ygrid); ss = xx+1i*yy; uss = u(ss);
levels = -log(2) + log(1.1:0.1:2);
hold off, contour(xgrid,ygrid,uss,levels,'k','linewidth',1.4)
ylim([-0.9,0.9]), grid on, axis equal, FS = 'fontsize';
title(['Equipotential curves for the Chebyshev ' ...
            'distribution = Bernstein ellipses'],FS,9)
print -dpdf 'equipot_chebyshev.pdf'
