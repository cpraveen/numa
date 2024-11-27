import graph;
size(300,300);

xaxis("$x$",Arrow);
yaxis("$y$",Arrow);

real x0 = 1.5;
pen pf = red+linewidth(2);
pen pt = blue+linewidth(1.5);

real f(real x) { return exp(0.5*x) - 1.5; }
real df(real x) { return 0.5*exp(0.5*x); }
real ft(real x) { return f(x0) + (x-x0)*df(x0); }

draw(graph(f,0.1,2.3,operator ..),pf);
draw(graph(ft,0.3,2.2,operator ..),pt);

draw((x0,0)--(x0,f(x0)),dashed);
label("$x_0$", (x0,0), S);
label("$y=f(x_0)+(x-x_0)f'(x_0)$", (2.0,1.0), E, pt);
label("$y=f(x)$", (2.3,f(2.3)), E, pf);

real x1 = x0 - f(x0)/df(x0);
label("$x_1$",(x1,0),SE);
