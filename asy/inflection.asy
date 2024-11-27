import graph;
size(350,0);

xaxis("$x$",Arrow);
yaxis("$f(x)$",Arrow);

real x0 = 4.05;
pen pf = red+linewidth(2);
pen pt = blue+linewidth(1.5);

real f(real x) { return tanh(x-3); }
real df(real x) { return 1.0 - tanh(x-3)^2; }
//real ft(real x) { return f(x0) + (x-x0)*df(x0); }

draw(graph(f,0.0,6.0,operator ..),pf);
//draw(graph(ft,0.3,2.2,operator ..),pt);

label("$x_0$", (x0,0), S);
//label("$f(x_0)+(x-x_0)f'(x_0)$", (2.0,1.0), E, pt);

draw((x0,0)--(x0,f(x0)),dashed);
real x1 = x0 - f(x0)/df(x0);
label("$x_1$",(x1,0),SE);
draw((x0,f(x0))--(x1,0), dashed);
draw((x1,0)--(x1,f(x1)),dashed);

real x2 = x1 - f(x1)/df(x1);
draw((x1,f(x1))--(x2,0),dashed);
label("$x_2$",(x2,0),S);
