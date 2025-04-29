#include "raylib.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

float cos_acos_3(float x) {
	x = sqrt(0.5+0.5*x);
	return x*(x*(x*(x*-0.008972+0.039071)-0.107074)+0.576975)+0.5;
}

double bezier_sdf(double px, double py, double Ax, double Ay, double Bx, double By, double Cx, double Cy) {
	double ax = Bx - Ax;
	double ay = By - Ay;
	double bx = Ax - 2.0*Bx + Cx;
	double by = Ay - 2.0*By + Cy;
	double cx = 2.0*ax;
	double cy = 2.0*ay;
	double dx = Ax - px;
	double dy = Ay - py;

	double kk = 1.0 / (bx * bx + by * by);
	double kx = kk * (ax * bx + ay * by);
	double ky = kk * (2.0 * (ax * ax + ay * ay) + (dx * bx + dy * by)) / 3.0;
	double kz = kk * (dx * ax + dy * ay);

	double res = 0.0;
	double sgn = 0.0;

	double p  = ky - kx * kx;
	double q  = kx * (2.0 * kx * kx - 3.0 * ky) + kz;
	double p3 = p * p * p;
	double q2 = q * q;
	double h  = q2 + 4.0 * p3;

	if (h >= 0.0) {
		h = sqrt(h);
		double xx = (h - q) / 2.0;
		double xy = (-h -q) / 2.0;

		#if 0
		// When p≈0 and p<0, h-q has catastrophic cancelation. So, we do
		// h=√(q²+4p³)=q·√(1+4p³/q²)=q·√(1+w) instead. Now we approximate
		// √ by a linear Taylor expansion into h≈q(1+½w) so that the q's
		// cancel each other in h-q. Expanding and simplifying further we
		// get x=vec2(p³/q,-p³/q-q). And using a second degree Taylor
		// expansion instead: x=vec2(k,-k-q) with k=(1-p³/q²)·p³/q
		if(abs(p) < 0.001) {
			float k = p3 / q;              // linear approx
			//float k = (1.0-p3/q2)*p3/q;  // quadratic approx
			x = vec2(k, -k - q);
		}
		#endif

		double xxsign = (xx < 0) ? -1 : 1;
		double xysign = (xy < 0) ? -1 : 1;

		double uvx = xxsign * cbrt(abs(xx));
		double uvy = xysign * cbrt(abs(xy));

		double t = uvx + uvy;

		// from NinjaKoala - single newton iteration to account for cancellation
		t -= (t*(t*t+3.0*p)+q)/(3.0*t*t+3.0*p);

		t = clamp(t - kx, 0.0, 1.0);
		double wx = dx + (cx + bx * t) * t;
		double wy = dy + (cy + by * t) * t;
		double outQx = wx + px;
		double outQy = wy + py;

		res = wx * wx + wy * wy;

		double tempx = (cx + 2.0 * bx * t);
		double tempy = (cy + 2.0 * by * t);
		sgn = tempx * wy - tempy * wx;
	}
	else
	{   // 3 roots
		double z = sqrt(-p);
		#if 0
		double v = acos(q/(p*z*2.0))/3.0;
		double m = cos(v);
		double n = sin(v);
		#else
		double m = cos_acos_3(q/(p*z*2.0));
		double n = sqrt(1.0-m*m);
		#endif
		n *= sqrt(3.0);

		double tx = clamp((m+m)*z-kx, 0.0, 1.0);
		double ty = clamp((-n-m)*z-kx, 0.0, 1.0);
		double tz = clamp((n-m)*z-kx, 0.0, 1.0);

		double qxx = dx+(cx+bx*tx)*tx;
		double qxy = dy+(cy+by*tx)*tx;
		double d2x = qxx * qxx + qxy * qxy;
		double tempx = (ax+bx*tx);
		double tempy = (ay+by*tx);
		double sx = tempx*qxy-tempy*qxx;

		double qyx = dx+(cx+bx*ty)*ty;
		double qyy = dy+(cy+by*ty)*ty;
		double d2y = qyx * qyx + qyy * qyy;
		double temp2x = (ax+bx*ty);
		double temp2y = (ay+by*ty);
		double sy = temp2x*qyy-temp2y*qyx;

		if (d2x < d2y) {
			res = d2x;
			sgn = sx;
			double outQx = qxx + px;
			double outQy = qxy + py;
		} else {
			res = d2y;
			sgn = sy;
			double outQx = qyx + px;
			double outQy = qyy + py;
		}
	}

	return sqrt(res) * ((sgn < 0) ? -1 : 1);
}

double bezier_deriv(double t, double A, double B, double C) {
	return 2 * (1 - t) * (B - A) + 2 * t * (C - B);
}

double bezier_deriv2(double t, double A, double B, double C) {
	return 2*(A - 2*B + C);
}

double distance_beziers(double Ax1, double Ay1, double Bx1, double By1, double Cx1, double Cy1, double Ax2, double Ay2, double Bx2, double By2, double Cx2, double Cy2) {
	// converting bezier to normal form: at^2 + bt + c
	// (A - 2B + C)t^2 + (2B - 2A)t + A
	// function we are trying to find root of is bezier_sdf((ax*t^2 + bx*t + cx), (ay*t^2 + by*t + cy), Ax2, Ay2, Bx2, By2, Cx2, Cy2)
	double ax = Ax1 - 2*Bx1 + Cx1;
	double ay = Ay1 - 2*By1 + Cy1;
	double bx = 2*Bx1 - 2*Ax1;
	double by = 2*By1 - 2*Ay1;
	double cx = Ax1;
	double cy = Ay1;

	// secant method
	double t0 = 0;
	double t1 = 0.5;
	double t = 0;
	for (int i = 0; i < 4; i++) {
		double f_t0 = bezier_sdf((ax*t0*t0 + bx*t0 + cx), (ay*t0*t0 + by*t0 + cy), Ax2, Ay2, Bx2, By2, Cx2, Cy2);
		double f_t1 = bezier_sdf((ax*t1*t1 + bx*t1 + cx), (ay*t1*t1 + by*t1 + cy), Ax2, Ay2, Bx2, By2, Cx2, Cy2);
		t = t1 - f_t1*((t1 - t0)/(f_t1 - f_t0));
		t0 = t1;
		t1 = t;
		// cout << t << endl;
	}
	return t;
}

int main(void) {
	// cout << bezier_sdf(1, 0, 1, 1, 2, 0, 3, 3) << endl;

	// Intitialize Window
	// SetConfigFlags(FLAG_MSAA_4X_HINT);
	// SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	SetConfigFlags(FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	InitWindow(800, 450, "Collider");

	SetTargetFPS(60000);
	while (!WindowShouldClose()) {
		BeginDrawing();

		ClearBackground(RAYWHITE);

		Vector2 points1[3];
		points1[0] = Vector2{0, 0};
		points1[1] = Vector2{800, 0};
		points1[2] = Vector2{800, 450};
		DrawSplineBezierQuadratic(points1, 3, 1, RED);
		Vector2 points2[3];
		points2[0] = Vector2{800, 0};
		points2[1] = Vector2{0, 0};
		points2[2] = Vector2{0, 450};
		DrawSplineBezierQuadratic(points2, 3, 1, RED);

		double t = distance_beziers(points1[0].x, points1[0].y, points1[1].x, points1[1].y, points1[2].x, points1[2].y, points2[0].x, points2[0].y, points2[1].x, points2[1].y, points2[2].x, points2[2].y);

		double ax = points1[0].x - 2*points1[1].x + points1[2].x;
		double ay = points1[0].y - 2*points1[1].y + points1[2].y;
		double bx = 2*points1[1].x - 2*points1[0].x;
		double by = 2*points1[1].y - 2*points1[0].y;
		double cx = points1[0].x;
		double cy = points1[0].y;
		DrawCircle(ax*t*t + bx*t + cx, ay*t*t + by*t + cy, 10, GREEN);

		DrawText(TextFormat("%f", bezier_sdf(GetMousePosition().x, GetMousePosition().y, 0, 0, 800, 0, 800, 450)), 25, GetScreenHeight() - 50, 20, LIGHTGRAY);
		DrawText(TextFormat("%f", t), 25, GetScreenHeight() - 25, 20, LIGHTGRAY);

		DrawFPS(GetScreenHeight() - 25, GetScreenHeight() - 25);

		EndDrawing();
	}

	CloseWindow();
	return 0;
}