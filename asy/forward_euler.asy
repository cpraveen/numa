// Forward Euler Method Illustration for y' = y with exact solution y = exp(t)

import graph;
import math;
import settings;

// Page settings
settings.outformat = "pdf";
size(200, 150, IgnoreAspect);

// --- Define the problem ---
// ODE: y' = f(t, y) = y
// Exact solution: y(t) = exp(t)

real f(real t, real y) {
    return y;
}

// Euler step function
pair eulerStep(real t, real y, real h) {
    return (t + h, y + h * f(t, y));
}

// --- Parameters ---
real t0 = 0.3;      // initial time
real y0 = 1;        // initial value y(t0) = 1
real h  = 1.0;      // step size (increased for visibility)

// --- Compute points ---
pair p0 = (t0, y0);
pair p1 = eulerStep(t0, y0, h);
pair p2 = eulerStep(p1.x, p1.y, h);

// --- Exact solutions ---
real exact(real t) { return y0 * exp(t-t0);}

// Exact solution passing through (t_1, y_1)
real exact1(real t) { return p1.y * exp(t - p1.x);}

// Exact solution passing through (t_2, y_2)
real exact2(real t) { return p2.y * exp(t - p2.x);}

// --- Axis ranges ---
real tmin = -0.3, tmax = 3.0;
real ymin = -0.2, ymax = 8.5;

// --- Draw axes ---
xaxis("$t$", arrow=Arrow);
yaxis("$y$", arrow=Arrow);

// --- Mark t_0, t_1, t_2 on x-axis ---
draw((t0, 0) -- (t0, -0.1), black);
label("$t_0$", (t0, -0.1), S);

draw((p1.x, 0) -- (p1.x, -0.1), black);
label("$t_1$", (p1.x, -0.1), S);

draw((p2.x, 0) -- (p2.x, -0.1), black);
label("$t_2$", (p2.x, -0.1), S);

// --- Mark y_0, y_1, y_2 on y-axis ---
draw((0, y0) -- (-0.1, y0), black);
label("$y_0$", (-0.1, y0), W);

draw((0, p1.y) -- (-0.1, p1.y), black);
label("$y_1$", (-0.1, p1.y), W);

draw((0, p2.y) -- (-0.1, p2.y), black);
label("$y_2$", (-0.1, p2.y), W);

// --- Draw exact solution through (t_0, y_0) ---
draw(graph(exact, t0, 2.3), blue, "Exact $y=e^t$");

// --- Draw exact solution through (t_1, y_1) ---
draw(graph(exact1, p1.x, 2.6), green, "Exact through $(t_1,y_1)$");

// --- Draw exact solution through (t_2, y_2) ---
draw(graph(exact2, p2.x, 2.9), purple, "Exact through $(t_2,y_2)$");

// --- Draw tangent line at (t0, y0) ---
real slope0 = f(t0, y0);
pair tanEnd = (p1.x + 0.3, p1.y + 0.3 * slope0);
draw(p0 -- tanEnd, orange+dashed);

// --- Draw tangent line at (t1, y1) ---
real slope1 = f(p1.x, p1.y);
pair tan1Start = (p1.x - 0.2, p1.y - 0.2 * slope1);
pair tan1End   = (p2.x + 0.2, p2.y + 0.2 * slope1);
draw(tan1Start -- tan1End, black+dashed);

// --- Draw Euler steps (horizontal + vertical) ---
draw(p0 -- (p1.x, p0.y), red, Arrow);
draw((p1.x, p0.y) -- p1, red, Arrow);
draw(p1 -- (p2.x, p1.y), red, Arrow);
draw((p2.x, p1.y) -- p2, red, Arrow);

// --- Label the step segments ---
label("$h$", ((p0.x + p1.x)/2, p0.y), S);
label("$h$", ((p1.x + p2.x)/2, p1.y), S);

// --- Mark points ---
dot(p0, black+4);
dot(p1, black+4);
dot(p2, black+4);

// --- Mark exact points for comparison ---
pair e1 = (p1.x, exact(p1.x));
pair e2 = (p2.x, exact(p2.x));
dot(e1, blue);
dot(e2, blue);

// --- Legend ---
//label("Forward Euler with $h=" + string(h) + "$", (1.5, 7.5), UnFill);
