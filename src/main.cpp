#include "raylib.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <complex>

using namespace std;

// === Translation of GLSL vector types. Additionally expanded to double type ===

struct vec2 {
	double x;
	double y;
	// Constructors
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
	// Constructors
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
struct bezier {
	vec2 A;
	vec2 B;
	vec2 C;
	bezier(vec2 A_in, vec2 B_in, vec2 C_in) : A(A_in), B(B_in), C(C_in) {}
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

// === Bezier SDF, credit Adam Simmons (https://www.shadertoy.com/view/ltXSDB) ===

// Test if point p crosses line (a, b), returns sign of result
double testCross(vec2 a, vec2 b, vec2 p) {
	return sign((b.y-a.y) * (p.x-a.x) - (b.x-a.x) * (p.y-a.y));
}

// Determine which side we're on (using barycentric parameterization)
double signBezier(vec2 p, vec2 A, vec2 B, vec2 C) {
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

// Solve cubic equation for roots (trigonometric solution for a depressed cubic)
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
double bezierSDF(vec2 p, bezier in) {
	vec2 A = in.A;
	vec2 B = in.B;
	vec2 C = in.C;

	// if abs(B * 2.0 - A - C) less than 1e-6, add 1e-4 to B
	B = mix(B + vec2(1e-4), B, step(1e-6, abs(B * 2.0 - A - C)));

	// Minimize f(t) = |p - b(t)|^2 (v^2 = dot(v, v))
	// => f(t) = |p|^2 - 2p \cdot b(t) + |b(t)|^2
	// => To minimize, find where f'(t) = 0
	// => Solve f'(t) = - 2p \cdot b'(t) + 2b(t) \cdot b'(t) = b'(t) \cdot 2 * (b(t) - p) = 0
	// => (b(t) - p) \cdot b'(t) = 0
	// (b(t) - p) \cdot b'(t) = 0
	// => ((A - 2B + C)t^2 + (2B - 2A)t + A - p) \cdot (2(A - 2B + C)t + (2B - 2A)) = 0
	vec2 a = B - A,
	b = A - B * 2.0 + C,
	c = a * 2.0,
	d = A - p;
	vec3 k = vec3(3.0*dot(a,b), 2.0*dot(a,a) + dot(d,b), dot(d,a)) / dot(b,b);
	vec3 t = clamp(solveCubic(k.x, k.y, k.z), 0.0, 1.0);
	vec2 pos = A + (c + b*t.x)*t.x;
	double dis = length(pos - p);
	pos = A + (c + b*t.y)*t.y;
	dis = min(dis, length(pos - p));
	pos = A + (c + b*t.z)*t.z;
	dis = min(dis, length(pos - p));
	return dis * signBezier(p, A, B, C);
}

// Find closest point on bezier curve to point p
vec2 bezierClosest(vec2 p, bezier in) {
	vec2 A = in.A;
	vec2 B = in.B;
	vec2 C = in.C;
	B = mix(B + vec2(1e-4), B, step(1e-6, abs(B * 2.0 - A - C)));
	vec2 a = B - A,
	b = A - B * 2.0 + C,
	c = a * 2.0,
	d = A - p;
	vec3 k = vec3(3.0*dot(a,b), 2.0*dot(a,a) + dot(d,b), dot(d,a)) / dot(b,b);
	vec3 t = clamp(solveCubic(k.x, k.y, k.z), 0.0, 1.0);
	vec2 pos = A + (c + b*t.x)*t.x;
	return pos;
}

vec2 findAnalyticIntersection(bezier b1, bezier b2, int method, int n, double estimate0, double estimate1, double estimate2) {
	// Analytically find point on bezier A that causes SDF of bezier B to equal 0
	//
	// Methods:
	// 0 = secant method
	// 1 = bisection method
	// 2 = muller's method
	//
	// n: number of iterations to perform. 0 adjusts until accurate to 0.1, or a limit of 100
	//
	// Returns point of collision as a vec2

	// converting bezier to normal form: at^2 + bt + c
	// (A - 2B + C)t^2 + (2B - 2A)t + A
	vec2 a = b1.A - b1.B*2.0 + b1.C;
	vec2 b = (b1.B - b1.A)*2.0;
	vec2 c = b1.A;

	if (n == 0) n = 50;

	switch (method) {
		case 0: {
			// Secant method
			double t0 = estimate0;
			double t1 = estimate1;
			double t = 0;
			for (int i = 0; i < n; i++) {
				double f_t1 = bezierSDF((a*t1*t1 + b*t1 + c), b2);
				if (abs(f_t1) < 0.1) break;
				double f_t0 = bezierSDF((a*t0*t0 + b*t0 + c), b2);
				t = t1 - f_t1*((t1 - t0)/(f_t1 - f_t0));
				t0 = t1;
				t1 = t;
			}
			return a*t*t + b*t + c;
		}
		case 1: {
			// Bisection method
			// We must find t0, t1 such that f(t0) is negative and f(t1) is positive.
			double t0 = estimate0;
			double t1 = 1;
			double t = 0;
			double step = 2;
			double f_t1 = bezierSDF((a*t1*t1 + b*t1 + c), b2);
			for (int i = 0; i < pow(2, 5); i++) {
				double f_t0 = bezierSDF((a*t0*t0 + b*t0 + c), b2);
				if (f_t0 * f_t1 <= 0) break;
				if (t0 + step >= 1.0) {
					step /= 2.0;
					t0 = step / 2.0;
				}
				else t0 += step;
			}
			// Actual method
			for (int i = 0; i < n; i++) {
				t = (t0 + t1)/2;
				if (bezierSDF((a*t*t + b*t + c), b2) * bezierSDF((a*t0*t0 + b*t0 + c), b2) > 0) t0 = t;
				else t1 = t;
			}
			return a*t*t + b*t + c;
		}
		case 2: {
			// Muller's method
			double t0 = estimate0;
			double t1 = estimate1;
			double t2 = estimate2;
			for (int i = 0; i < n; i++) {
				double h0 = t1 - t0;
				double h1 = t2 - t1;
				double f_t0 = bezierSDF((a*t0*t0 + b*t0 + c), b2);
				double f_t1 = bezierSDF((a*t1*t1 + b*t1 + c), b2);
				double f_t2 = bezierSDF((a*t2*t2 + b*t2 + c), b2);
				double d0 = (f_t1 - f_t0)/h0;
				double d1 = (f_t2 - f_t1)/h1;
				double a_m = (d1 - d0)/(h1 + h0);
				double b_m = a_m*h1 + d1;
				double c_m = f_t2;
				double t3 = t2 - (2*c_m)/max(
					abs(b_m + pow(pow(b_m, 2) - 4*a_m*c_m, 1/2)),
					abs(b_m - pow(pow(b_m, 2) - 4*a_m*c_m, 1/2))
				);
				t0 = t1;
				t1 = t2;
				t2 = t3;
			}
			return a*t2*t2 + b*t2 + c;
		}
		default: {
			cout << "Illegal Method" << endl;
			return vec2(0, 0);
		}
	}
}

vector<vec2> findAllAnalyticIntersection(bezier b1, bezier b2, int method, int n) {
	// Results only include those with error less than 0.1
	// Points of intersection are considered separate if separated by distance of at least 1
	vector<vec2> results;
	int sections = 10;
	for (int i = 0; i <= sections; i++) {
		vec2 result = findAnalyticIntersection(b1, b2, method, n, i*(1.0/sections), 0.5, 1.0);
		int isIn = 0;
		for (vec2 r : results) {
			if (abs(r.x - result.x) < 1 && abs(r.y - result.y) < 1) {
				isIn = 1;
				break;
			}
		}
		if (isIn == 0 && abs(bezierSDF(result, b1)) < 0.1 && abs(bezierSDF(result, b2)) < 0.1) {
			results.push_back(result);
		}
	}
	return results;
}

struct mcresult {
	vec2 p;
	int in;
	mcresult(vec2 p_in, int in_in) : p(p_in), in(in_in) {}
};
vector<mcresult> montecarlo(vector<vector<bezier>> shapes, int x, int y, int width, int height, int n) {
	vector<mcresult> output;
	srand(time(0));
	for (int i = 0; i < n; i++) {
		vec2 rpoint = vec2(x + (rand() % (width - x)), y + (rand() % (height - y)));
		int num_in = 0;
		for (vector<bezier> shape : shapes) {
			int in = 1;
			for (bezier b : shape) {
				if (bezierSDF(rpoint, b) > 0) {
					in = 0;
					break;
				}
			}
			if (in) num_in++;
		}
		output.push_back(mcresult(rpoint, num_in));
	}
	return output;
}

int main(void) {
	// Intitialize Window
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	// SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	SetConfigFlags(FLAG_WINDOW_HIGHDPI | FLAG_WINDOW_RESIZABLE);
	InitWindow(800, 600, "Collider");
	SetTargetFPS(60);

	vector<vector<bezier>> shapes;
	shapes.push_back({bezier(vec2(20, 20), vec2(GetScreenWidth() - 220, 20), vec2(GetScreenWidth() - 220, GetScreenHeight() - 20))});
	shapes.push_back({bezier(vec2(20, GetScreenHeight() - 40), vec2(20, 20), vec2(GetScreenWidth() - 220, 20))});

	vec2* adjusting = nullptr;

	int n = 0;
	double p[3] = {0.0, 0.5, 1.0};
	int changing_p = 0;
	int method = 0;
	int show_sdf = 0;

	vector<mcresult> montecarlo_results;
	int montecarlo_n = 2000;

	const float ctrlpt_size = 10;

	while (!WindowShouldClose()) {
		BeginDrawing();
		ClearBackground(RAYWHITE);

		// Handle all shapes
		for (int i = 0; i < shapes.size(); i++) {
			for (int i2 = 0; i2 < shapes[i].size(); i2++) {
				bezier b = shapes[i][i2];
				// Draw control points
				DrawRectangle(b.A.x - ctrlpt_size/2, b.A.y - ctrlpt_size/2, ctrlpt_size, ctrlpt_size, ORANGE);
				DrawRectangle(b.B.x - ctrlpt_size/2, b.B.y - ctrlpt_size/2, ctrlpt_size, ctrlpt_size, ORANGE);
				DrawRectangle(b.C.x - ctrlpt_size/2, b.C.y - ctrlpt_size/2, ctrlpt_size, ctrlpt_size, ORANGE);
				// Identify control point adjustment
				if (IsMouseButtonUp(MOUSE_BUTTON_LEFT)) adjusting = nullptr;
				else if (adjusting == nullptr && abs(b.A.x - GetMousePosition().x) < 5 && abs(b.A.y - GetMousePosition().y) < 5) adjusting = &shapes[i][i2].A;
				else if (adjusting == nullptr && abs(b.B.x - GetMousePosition().x) < 5 && abs(b.B.y - GetMousePosition().y) < 5) adjusting = &shapes[i][i2].B;
				else if (adjusting == nullptr && abs(b.C.x - GetMousePosition().x) < 5 && abs(b.C.y - GetMousePosition().y) < 5) adjusting = &shapes[i][i2].C;
				// Draw contact points with all other shapes
				for (int i3 = i + 1; i3 < shapes.size(); i3++) {
					for (bezier b2 : shapes[i3]) {
						// for (vec2 p : findAllAnalyticIntersection(b, b2, 1, n)) {
						// 	DrawCircle(p.x, p.y, 5, GREEN);
						// }
						vec2 pt = findAnalyticIntersection(b, b2, method, n, p[0], p[1], p[2]);
						DrawCircle(pt.x, pt.y, 5, GREEN);
					}
				}
				// Draw the bezier
				Vector2 points[3] = {
					Vector2{static_cast<float>(b.A.x), static_cast<float>(b.A.y)},
					Vector2{static_cast<float>(b.B.x), static_cast<float>(b.B.y)},
					Vector2{static_cast<float>(b.C.x), static_cast<float>(b.C.y)}
				};
				DrawSplineBezierQuadratic(points, 3, 1, RED);
				// Draw line from mouse pos to closest point
				if (show_sdf) {
					vec2 closest_to_mouse = bezierClosest(vec2(GetMouseX(), GetMouseY()), b);
					DrawLineV(Vector2{static_cast<float>(closest_to_mouse.x), static_cast<float>(closest_to_mouse.y)}, GetMousePosition(), LIGHTGRAY);
				}
			}
		}

		// Perform control point adjustment
		if (adjusting != nullptr) {
			adjusting->x = GetMousePosition().x;
			adjusting->y = GetMousePosition().y;
		}

		// Keyboard controls
		if (IsKeyPressed(KEY_ONE)) changing_p = 0;
		if (IsKeyPressed(KEY_TWO)) changing_p = 1;
		if (IsKeyPressed(KEY_THREE)) changing_p = 2;
		if (IsKeyPressed(KEY_UP)) n++;
		if (IsKeyPressed(KEY_DOWN) && n > 0) n--;
		if (IsKeyPressed(KEY_RIGHT) && p[changing_p] < 0.99) p[changing_p] += 0.1;
		if (IsKeyPressed(KEY_LEFT) && p[changing_p] > 0.01) p[changing_p] -= 0.1;
		if (IsKeyPressed(KEY_SPACE)) show_sdf = !show_sdf;

		// === Handle UI ===
		DrawRectangle(GetScreenWidth() - 200, 0, 200, GetScreenHeight(), LIGHTGRAY);
		// Secant Button
		Rectangle secant_button = Rectangle{static_cast<float>(GetScreenWidth() - 175), 25, 150, 45};
		if (method == 0) DrawRectangleRec(secant_button, DARKGRAY);
		else DrawRectangleRec(secant_button, GRAY);
		DrawText("Secant", GetScreenWidth() - 165, 35, 25, BLACK);
		if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), secant_button)) {
			method = 0;
		}
		// Bisection Button
		Rectangle bisection_button = Rectangle{static_cast<float>(GetScreenWidth() - 175), 95, 150, 45};
		if (method == 1) DrawRectangleRec(bisection_button, DARKGRAY);
		else DrawRectangleRec(bisection_button, GRAY);
		DrawText("Bisection", GetScreenWidth() - 165, 105, 25, BLACK);
		if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), bisection_button)) {
			method = 1;
		}
		// Muller Button
		Rectangle muller_button = Rectangle{static_cast<float>(GetScreenWidth() - 175), 165, 150, 45};
		if (method == 2) DrawRectangleRec(muller_button, DARKGRAY);
		else DrawRectangleRec(muller_button, GRAY);
		DrawText("Muller", GetScreenWidth() - 165, 175, 25, BLACK);
		if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), muller_button)) {
			method = 2;
		}
		// Monte Carlo button
		Rectangle montecarlo_button = Rectangle{static_cast<float>(GetScreenWidth() - 175), 235, 150, 45};
		if (montecarlo_results.size() != 0) DrawRectangleRec(montecarlo_button, DARKGRAY);
		else DrawRectangleRec(montecarlo_button, GRAY);
		DrawText("M. Carlo", GetScreenWidth() - 165, 245, 25, BLACK);
		if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), montecarlo_button)) {
			if (montecarlo_results.size() == 0) montecarlo_results = montecarlo(shapes, 0, 0, GetScreenWidth() - 200, GetScreenHeight(), montecarlo_n);
			else montecarlo_results.clear();
		}
		// Monte Carlo results
		if (montecarlo_results.size() > 0) {
			int two_or_more = 0;
			for (mcresult r : montecarlo_results) {
				if (r.in == 0) DrawCircleV(Vector2{static_cast<float>(r.p.x), static_cast<float>(r.p.y)}, 2, YELLOW);
				else if (r.in == 1) DrawCircleV(Vector2{static_cast<float>(r.p.x), static_cast<float>(r.p.y)}, 2, PINK);
				else if (r.in >= 2) {
					DrawCircleV(Vector2{static_cast<float>(r.p.x), static_cast<float>(r.p.y)}, 2, BLUE);
					two_or_more++;
				}
			}
			DrawRectangle((GetScreenWidth() - 200) / 2 - 200, GetScreenHeight() - 30, 400, 45, LIGHTGRAY);
			float area = ((float) two_or_more / (float) montecarlo_n) * ((GetScreenWidth() - 200) * GetScreenHeight());
			DrawText(TextFormat("2+ shape overlap area: %i", (int) round(area)), (GetScreenWidth() - 200) / 2 - 190, GetScreenHeight() - 25, 20, DARKGRAY);
		}
		// n
		DrawText(TextFormat("n: %i", n), GetScreenWidth() - 175, 305, 25, BLACK);
		// p
		DrawText(TextFormat("p0: %f", p[0]), GetScreenWidth() - 175, 355, 25, changing_p == 0 ? GRAY : BLACK);
		DrawText(TextFormat("p1: %f", p[1]), GetScreenWidth() - 175, 405, 25, changing_p == 1 ? GRAY : BLACK);
		DrawText(TextFormat("p2: %f", p[2]), GetScreenWidth() - 175, 455, 25, changing_p == 2 ? GRAY : BLACK);
		// Instructions
		DrawText("SPACE: SDF", GetScreenWidth() - 175, GetScreenHeight() - 100, 20, BLACK);
		DrawText("1,2,3: select p", GetScreenWidth() - 175, GetScreenHeight() - 75, 20, BLACK);
		DrawText("U/D: change n", GetScreenWidth() - 175, GetScreenHeight() - 50, 20, BLACK);
		DrawText("L/R: change p", GetScreenWidth() - 175, GetScreenHeight() - 25, 20, BLACK);

		// DrawFPS(10, GetScreenHeight() - 25);

		EndDrawing();
	}

	CloseWindow();
	return 0;
}