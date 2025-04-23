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

int main(void)
{
	cout << bezier_sdf(1, 0, 1, 1, 2, 0, 3, 3) << endl;

	// Intitialize Window
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	InitWindow(800, 450, "Collider");

	Vector2 ballPosition = { GetScreenWidth()/2.0f, GetScreenHeight()/2.0f };
	Vector2 ballSpeed = { 5.0f, 4.0f };
	int ballRadius = 20;

	bool pause = 0;
	int framesCounter = 0;

	// Main game loop
	SetTargetFPS(60);
	while (!WindowShouldClose())
	{
		if (IsKeyPressed(KEY_SPACE)) pause = !pause;

		if (!pause)
		{
			ballPosition.x += ballSpeed.x;
			ballPosition.y += ballSpeed.y;

			// Check walls collision for bouncing
			if ((ballPosition.x >= (GetScreenWidth() - ballRadius)) || (ballPosition.x <= ballRadius)) ballSpeed.x *= -1.0f;
			if ((ballPosition.y >= (GetScreenHeight() - ballRadius)) || (ballPosition.y <= ballRadius)) ballSpeed.y *= -1.0f;
		}
		else framesCounter++;

		// Draw
		BeginDrawing();

			ClearBackground(RAYWHITE);

			DrawCircleV(ballPosition, (float)ballRadius, MAROON);
			DrawText("PRESS SPACE to PAUSE BALL MOVEMENT", 10, GetScreenHeight() - 25, 20, LIGHTGRAY);

			// On pause, we draw a blinking message
			if (pause && ((framesCounter/30)%2)) DrawText("PAUSED", 350, 200, 30, GRAY);

			DrawFPS(10, 10);

		EndDrawing();
	}

	CloseWindow();
	return 0;
}