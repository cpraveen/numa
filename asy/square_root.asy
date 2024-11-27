import graph;
size(300,300);

xaxis("$x$",Arrow);
yaxis("$y$",Arrow);

real a  = 2.0;
real x0 = 0.7;
pen pf = red+linewidth(2);
pen pt = blue+linewidth(1.5);

real f(real x) { return x*x - a; }
real df(real x) { return 2.0*x; }
real ft(real x) { return f(x0) + (x-x0)*df(x0); }

draw(graph(f,-2,2,operator ..),pf);
draw(graph(ft,0.3,2.4,operator ..),pt);

draw((x0,0)--(x0,f(x0)),dashed);
label("$x_0$", (x0,0), N);
label("$y=f(x_0)+(x-x_0)f'(x_0)$", (2.0,1.0), E, pt);
label("$y=x^2-a$", (2.0,f(2.0)), E, pf);

real x1 = x0 - f(x0)/df(x0);
label("$x_1$",(x1,0),SE);
label("$-a$",(-sqrt(a),0),NE);
label("$+a$",(+sqrt(a),0),NW);
