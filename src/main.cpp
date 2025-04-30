#include "raylib.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <array>

using namespace std;

// === Translation of GLSL vector types. Additionally expanded to double type ===
struct vec2 {
	double x;
	double y;
	// constructors
	vec2(double x_in, double y_in) : x(x_in), y(y_in) {}
	vec2(double both_in) : x(both_in), y(both_in) {}
	// Vector and Vector operations
	vec2 operator-(const vec2& other) const {
		return vec2(x-other.x, y-other.y);
	}
	vec2 operator+(const vec2& other) const {
		return vec2(x+other.x, y+other.y);
	}
	vec2 operator*(const vec2& other) const {
		return vec2(x*other.x, y*other.y);
	}
	// Vector and number operations
	vec2 operator+(double num) const {
		return vec2(x+num, y+num);
	}
	vec2 operator-(double num) const {
		return vec2(x-num, y-num);
	}
	vec2 operator/(double num) const {
		return vec2(x/num, y/num);
	}
	vec2 operator*(double num) const {
		return vec2(x*num, y*num);
	}
};
struct vec3 {
	double x;
	double y;
	double z;
	// constructors
	vec3(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}
	vec3(double both_in) : x(both_in), y(both_in), z(both_in) {}
	// Vector and number operations
	vec3 operator+(double num) const {
		return vec3(x+num, y+num, z+num);
	}
	vec3 operator*(double num) const {
		return vec3(x*num, y*num, z*num);
	}
	vec3 operator/(double num) const {
		return vec3(x/num, y/num, z/num);
	}
};

// === Translation of GLSL Functions ===

// Return sign of number or vector
double sign(double x) {return x >= 0 ? 1.0 : -1.0;}
vec2 sign(vec2 v) {return vec2(sign(v.x), sign(v.y));}

// Absolute value for vectors
vec2 abs(vec2 v) {return vec2(abs(v.x), abs(v.y));}

// Power function for two vectors
vec2 pow(vec2 v1, vec2 v2) {return vec2(pow(v1.x, v2.x), pow(v1.y, v2.y));}

// Linearly interpolate between x and y using weight a
double mix(double x, double y, double a) {return x * (1.0 - a) + y * a;}
vec2 mix(vec2 x, vec2 y, vec2 a) {return vec2(mix(x.x, y.x, a.x), mix(x.y, y.y, a.y));}

// Generate a step function by comparing two values
double step(double edge, double x) {return x < edge ? 0.0 : 1.0;}
vec2 step(vec2 edge, vec2 x) {return vec2(step(edge.x, x.x), step(edge.y, x.y));}

// Dot product of two vectors
double dot(vec2 a, vec2 b) {return a.x*b.x + a.y*b.y;}

// Clamping scalar or vector between min and max
double clamp(double x, double minVal, double maxVal) {return min(max(x, minVal), maxVal);}
vec3 clamp(vec3 v, double minVal, double maxVal) {return vec3(clamp(v.x, minVal, maxVal), clamp(v.y, minVal, maxVal), clamp(v.z, minVal, maxVal));}

// Length of vector
double length(vec2 v) {return sqrt(pow(v.x, 2) + pow(v.y, 2));}

// === Bezier SDF ===

// Test if point p crosses line (a, b), returns sign of result
double testCross(vec2 a, vec2 b, vec2 p) {
	return sign((b.y-a.y) * (p.x-a.x) - (b.x-a.x) * (p.y-a.y));
}

// Determine which side we're on (using barycentric parameterization)
double signBezier(vec2 A, vec2 B, vec2 C, vec2 p) {
	vec2 a = C - A, b = B - A, c = p - A;
	vec2 bary = vec2(c.x*b.y-b.x*c.y, a.x*c.y-c.x*a.y) / (a.x*b.y-b.x*a.y);
	vec2 d = vec2(bary.y * 0.5, 0.0) + 1.0 - bary.x - bary.y;
	return mix(
		sign(d.x * d.x - d.y),
		mix(-1.0, 1.0,
			step(
				testCross(A, B, p) * testCross(B, C, p),
				0.0
			)
		),
		step((d.x - d.y), 0.0)
	) * testCross(A, C, B);
}

// Solve cubic equation for roots
vec3 solveCubic(double a, double b, double c) {
	double p = b - a*a / 3.0, p3 = p*p*p;
	double q = a * (2.0*a*a - 9.0*b) / 27.0 + c;
	double d = q*q + 4.0*p3 / 27.0;
	double offset = -a / 3.0;
	if(d >= 0.0) {
		double z = sqrt(d);
		vec2 x = (vec2(z, -z) - q) / 2.0;
		vec2 uv = sign(x)*pow(abs(x), vec2(1.0/3.0));
		return vec3(offset + uv.x + uv.y);
	}
	double v = acos(-sqrt(-27.0 / p3) * q / 2.0) / 3.0;
	double m = cos(v), n = sin(v)*1.732050808;
	return vec3(m + m, -n - m, n - m) * sqrt(-p / 3.0) + offset;
}

// Find the signed distance from a point to a bezier curve
double sdBezier(vec2 p, vec2 A, vec2 B, vec2 C) {
	B = mix(B + vec2(1e-4), B, step(1e-6, abs(B * 2.0 - A - C)));
	vec2 a = B - A, b = A - B * 2.0 + C, c = a * 2.0, d = A - p;
	vec3 k = vec3(3.*dot(a,b),2.*dot(a,a)+dot(d,b),dot(d,a)) / dot(b,b);
	vec3 t = clamp(solveCubic(k.x, k.y, k.z), 0.0, 1.0);
	vec2 pos = A + (c + b*t.x)*t.x;
	double dis = length(pos - p);
	pos = A + (c + b*t.y)*t.y;
	dis = min(dis, length(pos - p));
	pos = A + (c + b*t.z)*t.z;
	dis = min(dis, length(pos - p));
	return dis * signBezier(A, B, C, p);
}

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

array<double, 2> collision_beziers(double Ax1, double Ay1, double Bx1, double By1, double Cx1, double Cy1, double Ax2, double Ay2, double Bx2, double By2, double Cx2, double Cy2, int method, int n, double starting_estimate) {
	// Ax1: x coordinate of P0 in bezier 1, Cy2: y coordinate of P2 in bezier 2, etc.
	//
	// Methods:
	// 0 = secant method
	//
	// n: number of iterations to perform. 0 adjusts until accurate, or a limit of 100
	//
	// Returns point of collision as an array<double,2>{x, y}

	// converting bezier to normal form: at^2 + bt + c
	// (A - 2B + C)t^2 + (2B - 2A)t + A
	// function we are trying to find root of is bezier_sdf((ax*t^2 + bx*t + cx), (ay*t^2 + by*t + cy), Ax2, Ay2, Bx2, By2, Cx2, Cy2)
	double ax = Ax1 - 2*Bx1 + Cx1;
	double ay = Ay1 - 2*By1 + Cy1;
	double bx = 2*Bx1 - 2*Ax1;
	double by = 2*By1 - 2*Ay1;
	double cx = Ax1;
	double cy = Ay1;

	if (n == 0) n = 100;

	// secant method
	double t0 = 0;
	double t1 = starting_estimate;
	double t = 0;
	for (int i = 0; i < n; i++) {
		double f_t1 = bezier_sdf((ax*t1*t1 + bx*t1 + cx), (ay*t1*t1 + by*t1 + cy), Ax2, Ay2, Bx2, By2, Cx2, Cy2);
		if (abs(f_t1) <= 0.1) break;
		double f_t0 = bezier_sdf((ax*t0*t0 + bx*t0 + cx), (ay*t0*t0 + by*t0 + cy), Ax2, Ay2, Bx2, By2, Cx2, Cy2);
		t = t1 - f_t1*((t1 - t0)/(f_t1 - f_t0));
		t0 = t1;
		t1 = t;
		// cout << t << endl;
	}
	return array<double, 2>{ax*t*t + bx*t + cx, ay*t*t + by*t + cy};
}

int main(void) {
	// Intitialize Window
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	// SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	SetConfigFlags(FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	InitWindow(1280, 720, "Collider");
	SetTargetFPS(60);

	Vector2 points1[3];
	points1[0] = Vector2{0, 0};
	points1[1] = Vector2{800, 0};
	points1[2] = Vector2{800, 450};
	Vector2 points2[3];
	points2[0] = Vector2{800, 0};
	points2[1] = Vector2{0, 0};
	points2[2] = Vector2{0, 450};

	Vector2* adjusting = nullptr;

	int n = 3;
	double starting_estimate = 0.5;

	while (!WindowShouldClose()) {
		BeginDrawing();
		ClearBackground(RAYWHITE);

		for (int i = 0; i < 3; i++) {
			DrawRectangle(points1[i].x - 5, points1[i].y - 5, 10, 10, ORANGE);
			if (IsMouseButtonUp(MOUSE_BUTTON_LEFT)) adjusting = nullptr;
			else if (adjusting == nullptr && abs(points1[i].x - GetMousePosition().x) < 5 && abs(points1[i].y - GetMousePosition().y) < 5) adjusting = &points1[i];
		}
		for (int i = 0; i < 3; i++) {
			DrawRectangle(points2[i].x - 5, points2[i].y - 5, 10, 10, ORANGE);
			if (IsMouseButtonUp(MOUSE_BUTTON_LEFT)) adjusting = nullptr;
			else if (adjusting == nullptr && abs(points2[i].x - GetMousePosition().x) < 5 && abs(points2[i].y - GetMousePosition().y) < 5) adjusting = &points2[i];
		}
		if (adjusting != nullptr) {
			(*adjusting).x = GetMousePosition().x;
			(*adjusting).y = GetMousePosition().y;
		}
		if (IsKeyPressed(KEY_UP)) n++;
		if (IsKeyPressed(KEY_DOWN) && n > 0) n--;
		if (IsKeyPressed(KEY_RIGHT) && starting_estimate < 0.99) starting_estimate += 0.1;
		if (IsKeyPressed(KEY_LEFT) && starting_estimate > 0.01) starting_estimate -= 0.1;

		DrawSplineBezierQuadratic(points1, 3, 1, RED);
		DrawSplineBezierQuadratic(points2, 3, 1, BLUE);

		array<double, 2> collision_point = collision_beziers(
			points1[0].x, points1[0].y, points1[1].x, points1[1].y, points1[2].x, points1[2].y,
			points2[0].x, points2[0].y, points2[1].x, points2[1].y, points2[2].x, points2[2].y,
			0, n, starting_estimate);

		double ax = points1[0].x - 2*points1[1].x + points1[2].x;
		double ay = points1[0].y - 2*points1[1].y + points1[2].y;
		double bx = 2*points1[1].x - 2*points1[0].x;
		double by = 2*points1[1].y - 2*points1[0].y;
		double cx = points1[0].x;
		double cy = points1[0].y;
		DrawCircle(collision_point[0], collision_point[1], 5, GREEN);

		// === Draw UI ===
		DrawRectangle(GetScreenWidth() - 200, 0, 200, GetScreenHeight(), LIGHTGRAY);
		// Secant Button
		DrawRectangle(GetScreenWidth() - 175, 25, 150, 45, GRAY);
		DrawText("Secant", GetScreenWidth() - 165, 35, 25, BLACK);
		// Bisection Button
		DrawRectangle(GetScreenWidth() - 175, 95, 150, 45, GRAY);
		DrawText("Bisection", GetScreenWidth() - 165, 105, 25, BLACK);
		// Muller Button
		DrawRectangle(GetScreenWidth() - 175, 165, 150, 45, GRAY);
		DrawText("Muller", GetScreenWidth() - 165, 175, 25, BLACK);
		// n
		DrawText(TextFormat("n: %i", n), GetScreenWidth() - 175, 235, 25, BLACK);
		// p0
		DrawText(TextFormat("p0: %f", starting_estimate), GetScreenWidth() - 175, 285, 25, BLACK);
		// Instructions
		DrawText("U/D: change n", GetScreenWidth() - 190, GetScreenHeight() - 50, 20, BLACK);
		DrawText("L/R: change p0", GetScreenWidth() - 190, GetScreenHeight() - 25, 20, BLACK);

		// DrawText(TextFormat("%f", bezier_sdf(GetMousePosition().x, GetMousePosition().y, points1[0].x, points1[0].y, points1[1].x, points1[1].y, points1[2].x, points1[2].y)), 25, GetScreenHeight() - 75, 20, LIGHTGRAY);
		DrawText(TextFormat("%f", sdBezier(vec2(GetMousePosition().x, GetMousePosition().y), vec2(points1[0].x, points1[0].y), vec2(points1[1].x, points1[1].y), vec2(points1[2].x, points1[2].y))), 25, GetScreenHeight() - 75, 20, LIGHTGRAY);
		DrawText(TextFormat("%f %f", collision_point[0], collision_point[1]), 25, GetScreenHeight() - 50, 20, LIGHTGRAY);
		DrawText(TextFormat("%i", n), 25, GetScreenHeight() - 25, 20, LIGHTGRAY);

		DrawFPS(GetScreenHeight() - 25, GetScreenHeight() - 25);

		EndDrawing();
	}

	CloseWindow();
	return 0;
}