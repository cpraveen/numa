import graph;
size(300,300);

xaxis("$x$",Arrow);
yaxis("$y$",Arrow);

real a = 0.5;
real x0 = 0.6;
pen pf = red+linewidth(2);
pen pt = blue+linewidth(1.5);

real f(real x) { return 1.0/x - a; }
real df(real x) { return -1.0/(x*x); }
real ft(real x) { return f(x0) + (x-x0)*df(x0); }

draw(graph(f,0.3,2.8,operator ..),pf);
draw(graph(ft,0.3,1.4,operator ..),pt);

draw((x0,0)--(x0,f(x0)),dashed);
label("$x_0$", (x0,0), S);
label("$y=f(x_0)+(x-x_0)f'(x_0)$", (1.3,ft(1.3)), E, pt);
label("$y=\frac{1}{x}-a$", (0.3,f(0.3)), E, pf);

real x1 = x0 - f(x0)/df(x0);
label("$x_1$",(x1,0),SW);

label("$\frac{1}{a}$",(1.0/a,0.0),S);
// asymptote
draw((0,-a)--(2.8,-a),dashed);
label("$-a$",(0,-a),W);
label("asymptote", (0.5,-a),N);
