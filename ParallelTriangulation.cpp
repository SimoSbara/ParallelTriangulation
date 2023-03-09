// ParallelTriangulation.cpp : Questo file contiene la funzione 'main', in cui inizia e termina l'esecuzione del programma.
//

#include <iostream>
#include <CDT.h>
#include <omp.h>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <Windows.h>
#include <chrono>

#include "imgui.h"
#include "imgui-SFML.h"

#pragma comment(lib, "sfml-system.lib")
#pragma comment(lib, "sfml-graphics.lib")
#pragma comment(lib, "sfml-window.lib")
#pragma comment(lib, "opengl32.lib")

#define NUM_THREADS 12
#define EXAMPLE_FOLDER "examples/"
#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600

struct Quadrante
{
	double xMin;
	double xMax;
	double yMin;
	double yMax;

	bool isInside(double x, double y)
	{
		return (xMin <= x && x <= xMax) && (yMin <= y && y <= yMax);
	}
};

CDT::Triangulation<double>* triangulator; //singolo thread
std::vector<CDT::V2d<double>> points;
std::vector<double> zHeight;

std::vector<double> zHeights[NUM_THREADS];
std::vector<CDT::V2d<double>> tPoints[NUM_THREADS];
CDT::Triangulation<double>* triangulators[NUM_THREADS];
Quadrante quadranti[NUM_THREADS];

std::string inputName;

sf::VertexArray trianglesToDraw(sf::Triangles);

double xMin = DBL_MAX;
double xMax = DBL_MIN;
double yMin = DBL_MAX;
double yMax = DBL_MIN;

HANDLE  threadTriangulation;
HANDLE  threadActions;
DWORD   threadTrID;
DWORD   threadActID;

bool triangulatioDone;
bool drawTriangles;
bool isRunning;

float zoom = 1.0;
float shiftX = 0;
float shiftY = 0;

typedef struct {
	double r;       // a fraction between 0 and 1
	double g;       // a fraction between 0 and 1
	double b;       // a fraction between 0 and 1
} rgb;

typedef struct {
	double h;       // angle in degrees
	double s;       // a fraction between 0 and 1
	double v;       // a fraction between 0 and 1
} hsv;

//da tinta a rgb
rgb hsv2rgb(hsv in)
{
	double      hh, p, q, t, ff;
	long        i;
	rgb         out;

	if (in.s <= 0.0) {       // < is bogus, just shuts up warnings
		out.r = in.v;
		out.g = in.v;
		out.b = in.v;
		return out;
	}
	hh = in.h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = in.v * (1.0 - in.s);
	q = in.v * (1.0 - (in.s * ff));
	t = in.v * (1.0 - (in.s * (1.0 - ff)));

	switch (i) {
	case 0:
		out.r = in.v;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = in.v;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = in.v;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = in.v;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = in.v;
		break;
	case 5:
	default:
		out.r = in.v;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}

rgb HSVTable[360];

void DrawTriangles(double minXGlobal, double maxXGlobal, double minYGlobal, double maxYGlobal)
{
	std::vector<CDT::V2d<double>>& curPoints = points;
	double xDiff = 1.0 / (maxXGlobal - minXGlobal);
	double yDiff = 1.0 / (maxYGlobal - minYGlobal);

	//double curHSV = i * hsvThread;

	for (int j = 0; j < triangulator->triangles.size(); j++)
	{
		sf::Vertex v1;
		sf::Vertex v2;
		sf::Vertex v3;

		//sf::VertexArray triangleSFML(sf::LinesStrip, 3);
		sf::Vertex lineA[2];
		sf::Vertex lineB[2];
		sf::Vertex lineC[2];

		CDT::Triangle* triangle = &triangulator->triangles[j];

		//lineA[0] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[0]].x, curPoints[triangle->vertices[0]].y)); //ordinato antiorario
		//lineA[1] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[1]].x, curPoints[triangle->vertices[1]].y)); //ordinato antiorario

		//lineB[0] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[1]].x, curPoints[triangle->vertices[1]].y)); //ordinato antiorario
		//lineB[1] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[2]].x, curPoints[triangle->vertices[2]].y)); //ordinato antiorario

		//lineC[0] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[2]].x, curPoints[triangle->vertices[2]].y)); //ordinato antiorario
		//lineC[1] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[0]].x, curPoints[triangle->vertices[0]].y)); //ordinato antiorario


		// define the position of the triangle's points

		double normalizedX1 = (curPoints[triangle->vertices[0]].x - minXGlobal) * xDiff;
		double normalizedY1 = (curPoints[triangle->vertices[0]].y - minYGlobal) * yDiff;

		double normalizedX2 = (curPoints[triangle->vertices[1]].x - minXGlobal) * xDiff;
		double normalizedY2 = (curPoints[triangle->vertices[1]].y - minYGlobal) * yDiff;

		double normalizedX3 = (curPoints[triangle->vertices[2]].x - minXGlobal) * xDiff;
		double normalizedY3 = (curPoints[triangle->vertices[2]].y - minYGlobal) * yDiff;

		double X1, X2, X3;
		double Y1, Y2, Y3;
		double HSV1, HSV2, HSV3;

		double zoom = 1;
		double shiftX = 0;
		double shiftY = 0;

		shiftX *= zoom;
		shiftY *= zoom;

		X1 = normalizedX1 * (double)WINDOW_WIDTH;
		X2 = normalizedX2 * (double)WINDOW_WIDTH;
		X3 = normalizedX3 * (double)WINDOW_WIDTH;

		Y1 = normalizedY1 * (double)WINDOW_HEIGHT;
		Y2 = normalizedY2 * (double)WINDOW_HEIGHT;
		Y3 = normalizedY2 * (double)WINDOW_HEIGHT;

		HSV1 = normalizedX1 * 360.0;
		HSV2 = normalizedX2 * 360.0;
		HSV3 = normalizedX3 * 360.0;

		X1 *= zoom;
		X2 *= zoom;
		X3 *= zoom;

		Y1 *= zoom;
		Y2 *= zoom;
		Y3 *= zoom;

		X1 -= shiftX;
		X2 -= shiftX;
		X3 -= shiftX;

		Y1 -= shiftY;
		Y2 -= shiftY;
		Y3 -= shiftY;

		v1.position = sf::Vector2f(X1, Y1);
		v2.position = sf::Vector2f(X2, Y2);
		v3.position = sf::Vector2f(X3, Y3);

		// define the color of the triangle's points

		sf::Uint8 r1 = HSVTable[(int)HSV1].r * 255;
		sf::Uint8 g1 = HSVTable[(int)HSV1].g * 255;
		sf::Uint8 b1 = HSVTable[(int)HSV1].b * 255;

		sf::Uint8 r2 = HSVTable[(int)HSV2].r * 255;
		sf::Uint8 g2 = HSVTable[(int)HSV2].g * 255;
		sf::Uint8 b2 = HSVTable[(int)HSV2].b * 255;

		sf::Uint8 r3 = HSVTable[(int)HSV3].r * 255;
		sf::Uint8 g3 = HSVTable[(int)HSV3].g * 255;
		sf::Uint8 b3 = HSVTable[(int)HSV3].b * 255;


		v1.color = sf::Color(r1, g1, b1);
		v2.color = sf::Color(r2, g2, b2);
		v3.color = sf::Color(r3, g3, b3);

		v1.color = sf::Color::Red;
		v2.color = sf::Color::Green;
		v3.color = sf::Color::Blue;

		lineA[0].position = sf::Vector2f(X1, Y1);
		lineA[0].color = sf::Color::White;
		lineA[1].position = sf::Vector2f(X2, Y2);
		lineA[1].color = sf::Color::White;

		lineB[0].position = sf::Vector2f(X2, Y2);
		lineB[0].color = sf::Color::White;
		lineB[1].position = sf::Vector2f(X3, Y3);
		lineB[1].color = sf::Color::White;

		lineC[0].position = sf::Vector2f(X3, Y3);
		lineC[0].color = sf::Color::White;
		lineC[1].position = sf::Vector2f(X1, Y1);
		lineC[1].color = sf::Color::White;

		//window->draw(lineA, 2, sf::Lines);
		//window->draw(lineB, 2, sf::Lines);
		//window->draw(lineC, 2, sf::Lines);

		//window->draw(triangleSFML);
		//trianglesToDraw.push_back(triangleSFML);

		trianglesToDraw.append(v1);
		trianglesToDraw.append(v2);
		trianglesToDraw.append(v3);
	}

	triangulatioDone = TRUE;
}

void DrawTrianglesThreads(sf::RenderWindow* window, double minXGlobal, double maxXGlobal, double minYGlobal, double maxYGlobal)
{
	double xDiff = 1.0 / (maxXGlobal - minXGlobal);
	double yDiff = 1.0 / (maxYGlobal - minYGlobal);

	double hsvThread = 360.0 / (double)NUM_THREADS;

	for (int i = 0; i < NUM_THREADS; i++)
	{
		std::vector<CDT::V2d<double>>& curPoints = tPoints[i];

		double curHSV = i * hsvThread;

		for (int j = 0; j < triangulators[i]->triangles.size(); j++)
		{
			sf::VertexArray triangleSFML (sf::Triangles, 3);

			//sf::VertexArray triangleSFML(sf::LinesStrip, 3);
			sf::Vertex lineA[2];
			sf::Vertex lineB[2];
			sf::Vertex lineC[2];

			CDT::Triangle* triangle = &triangulators[i]->triangles[j];

			//lineA[0] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[0]].x, curPoints[triangle->vertices[0]].y)); //ordinato antiorario
			//lineA[1] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[1]].x, curPoints[triangle->vertices[1]].y)); //ordinato antiorario

			//lineB[0] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[1]].x, curPoints[triangle->vertices[1]].y)); //ordinato antiorario
			//lineB[1] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[2]].x, curPoints[triangle->vertices[2]].y)); //ordinato antiorario

			//lineC[0] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[2]].x, curPoints[triangle->vertices[2]].y)); //ordinato antiorario
			//lineC[1] = sf::Vertex(sf::Vector2f(curPoints[triangle->vertices[0]].x, curPoints[triangle->vertices[0]].y)); //ordinato antiorario


			// define the position of the triangle's points

			double normalizedX1 = (curPoints[triangle->vertices[0]].x - minXGlobal) * xDiff;
			double normalizedY1 = (curPoints[triangle->vertices[0]].y - minYGlobal) * yDiff;

			double normalizedX2 = (curPoints[triangle->vertices[1]].x - minXGlobal) * xDiff;
			double normalizedY2 = (curPoints[triangle->vertices[1]].y - minYGlobal) * yDiff;

			double normalizedX3 = (curPoints[triangle->vertices[2]].x - minXGlobal) * xDiff;
			double normalizedY3 = (curPoints[triangle->vertices[2]].y - minYGlobal) * yDiff;

			double X1, X2, X3;
			double Y1, Y2, Y3;
			double HSV1, HSV2, HSV3;

			double zoom = 1.0;
			double shiftX = 0;
			double shiftY = 0;

			shiftX *= zoom;
			shiftY *= zoom;

			X1 = normalizedX1 * 800.0;
			X2 = normalizedX2 * 800.0;
			X3 = normalizedX3 * 800.0;

			Y1 = normalizedY1 * 600.0;
			Y2 = normalizedY2 * 600.0;
			Y3 = normalizedY2 * 600.0;

			HSV1 = normalizedX1 * 360.0;
			HSV2 = normalizedX2 * 360.0;
			HSV3 = normalizedX3 * 360.0;

			X1 *= zoom;
			X2 *= zoom;
			X3 *= zoom;

			Y1 *= zoom;
			Y2 *= zoom;
			Y3 *= zoom;

			X1 -= shiftX;
			X2 -= shiftX;
			X3 -= shiftX;

			Y1 -= shiftY;
			Y2 -= shiftY;
			Y3 -= shiftY;

			triangleSFML[0].position = sf::Vector2f(X1, Y1);
			triangleSFML[1].position = sf::Vector2f(X2, Y2);
			triangleSFML[2].position = sf::Vector2f(X3, Y3);

			// define the color of the triangle's points

			sf::Uint8 r1 = HSVTable[(int)HSV1].r * 255;
			sf::Uint8 g1 = HSVTable[(int)HSV1].g * 255;
			sf::Uint8 b1 = HSVTable[(int)HSV1].b * 255;

			sf::Uint8 r2 = HSVTable[(int)HSV2].r * 255;
			sf::Uint8 g2 = HSVTable[(int)HSV2].g * 255;
			sf::Uint8 b2 = HSVTable[(int)HSV2].b * 255;

			sf::Uint8 r3 = HSVTable[(int)HSV3].r * 255;
			sf::Uint8 g3 = HSVTable[(int)HSV3].g * 255;
			sf::Uint8 b3 = HSVTable[(int)HSV3].b * 255;


			triangleSFML[0].color = sf::Color(r1, g1, b1);
			triangleSFML[1].color = sf::Color(r2, g2, b2);
			triangleSFML[2].color = sf::Color(r3, g3, b3);			
			
			triangleSFML[0].color = sf::Color::Red;
			triangleSFML[1].color = sf::Color::Green;
			triangleSFML[2].color = sf::Color::Blue;

			lineA[0].position = sf::Vector2f(X1, Y1);
			lineA[0].color = sf::Color::Red;
			lineA[1].position = sf::Vector2f(X2, Y2);
			lineA[1].color = sf::Color::Red;

			lineB[0].position = sf::Vector2f(X2, Y2);
			lineB[0].color = sf::Color::Red;
			lineB[1].position = sf::Vector2f(X3, Y3);
			lineB[1].color = sf::Color::Red;

			lineC[0].position = sf::Vector2f(X3, Y3);
			lineC[0].color = sf::Color::Red;
			lineC[1].position = sf::Vector2f(X1, Y1);
			lineC[1].color = sf::Color::Red;

			//window->draw(lineA, 2, sf::Lines);
			//window->draw(lineB, 2, sf::Lines);
			//window->draw(lineC, 2, sf::Lines);

			window->draw(triangleSFML);
		}
	}
}

void createOBJ()
{
	FILE* fp;
	fopen_s(&fp, "mesh.obj", "wb");

	if (fp)
	{
		for (int i = 0; i < triangulator->vertices.size(); i++)
		{
			fprintf(fp, "v %f %f %f\n", triangulator->vertices[i].x, triangulator->vertices[i].y, zHeight[i]);
		}

		for (int i = 0; i < triangulator->triangles.size(); i++)
		{
			fprintf(fp, "f %d %d %d\n", 
				triangulator->triangles[i].vertices[0] + 1, 
				triangulator->triangles[i].vertices[1] + 1,
				triangulator->triangles[i].vertices[2] + 1);
		}

		fclose(fp);
	}
}

DWORD WINAPI computeTriangulation(LPVOID lpParam)
{
	int npoints = 0;
	double x, y, z;

	FILE* fp;
	//fopen_s(&fp, "output.xyz", "rb");
	//fopen_s(&fp, "2023_03_05-23-22-21.xyz", "rb");
	fopen_s(&fp, inputName.c_str(), "rb");

	if (fp)
	{
		printf("leggo da file\n");

		fscanf_s(fp, "%d\n", &npoints);

		for (int i = 0; i < npoints; i++)
		{
			fscanf_s(fp, "%lf %lf %lf \n", &x, &y, &z);

			points.push_back({ x, y });
			zHeight.push_back(z);

			if (xMin > x)
				xMin = x;
			else if (xMax < x)
				xMax = x;

			if (yMin > y)
				yMin = y;
			else if (yMax < y)
				yMax = y;
		}

		fclose(fp);
	}

	{
		printf("separo i punti\n");

		double latoQuadrante = abs(xMax - xMin) / (double)NUM_THREADS;

		for (int i = 0; i < NUM_THREADS; i++)
		{
			double start = xMin + latoQuadrante * i;
			double end = start + latoQuadrante;

			quadranti[i].xMin = start;
			quadranti[i].xMax = end;

			quadranti[i].yMin = yMin;
			quadranti[i].yMax = yMax;
		}

		CDT::V2d<double> curPoint;

		for (int i = 0; i < npoints; i++)
		{
			curPoint = points[i];

			for (int j = 0; j < NUM_THREADS; j++)
			{
				if (quadranti[j].isInside(curPoint.x, curPoint.y))
				{
					tPoints[j].push_back(curPoint);
					zHeights[j].push_back(zHeight[i]);
					break;
				}
			}
		}

		printf("inizio single thread\n");
		auto start = omp_get_wtime();
		{
			triangulator = new CDT::Triangulation<double>(
				CDT::VertexInsertionOrder::Auto,
				CDT::IntersectingConstraintEdges::Resolve,
				1e-8);
			triangulator->insertVertices(points);
			triangulator->eraseSuperTriangle();
			//triangulator->eraseOuterTrianglesAndHoles();
			//CDT::RemoveDuplicates(triangulator->vertices);
		}
		auto end = omp_get_wtime();

		std::cout << "sequenziale tempo: " << end - start << " sec" << std::endl;

		createOBJ();

		printf("inizio parallelo\n");

		start = omp_get_wtime();

#pragma omp parallel num_threads(12) shared(triangulators, tPoints)
		{
			int curThread = omp_get_thread_num();

			triangulators[curThread] = new CDT::Triangulation<double>(
				CDT::VertexInsertionOrder::Auto,
				CDT::IntersectingConstraintEdges::Resolve,
				1e-9);

			triangulators[curThread]->insertVertices(tPoints[curThread]);
			triangulators[curThread]->eraseSuperTriangle();
			//triangulators[curThread]->eraseOuterTrianglesAndHoles();
			//CDT::RemoveDuplicates(triangulators[curThread]->vertices);
		}

		end = omp_get_wtime();

		std::cout << "parallelo: " << end - start << " sec" << std::endl;
	}

	drawTriangles = TRUE;
	while (drawTriangles)
	{
		Sleep(1);
	}

	triangulatioDone = TRUE;

	return 0;
}

DWORD WINAPI actionThread(LPVOID lpParam)
{
	std::vector<CDT::V2d<double>>& curPoints = points;

	//double curHSV = i * hsvThread;

	while (isRunning)
	{
		if (drawTriangles)
		{
			double xDiff = 1.0 / (xMax - xMin);
			double yDiff = 1.0 / (yMax - yMin);
			trianglesToDraw.clear();
			double actShiftX = shiftX * zoom;
			double actShiftY = shiftY * zoom;

			int numTriangles = triangulator->triangles.size();

#pragma omp parallel for
			for (int j = 0; j < numTriangles; j++)
			{
				sf::Vertex v1;
				sf::Vertex v2;
				sf::Vertex v3;

				//sf::VertexArray triangleSFML(sf::LinesStrip, 3);
				sf::Vertex lineA[2];
				sf::Vertex lineB[2];
				sf::Vertex lineC[2];

				CDT::Triangle* triangle = &triangulator->triangles[j];

				// define the position of the triangle's points

				double normalizedX1 = (curPoints[triangle->vertices[0]].x - xMin) * xDiff;
				double normalizedY1 = (curPoints[triangle->vertices[0]].y - yMin) * yDiff;

				double normalizedX2 = (curPoints[triangle->vertices[1]].x - xMin) * xDiff;
				double normalizedY2 = (curPoints[triangle->vertices[1]].y - yMin) * yDiff;

				double normalizedX3 = (curPoints[triangle->vertices[2]].x - xMin) * xDiff;
				double normalizedY3 = (curPoints[triangle->vertices[2]].y - yMin) * yDiff;

				double X1, X2, X3;
				double Y1, Y2, Y3;
				double HSV1, HSV2, HSV3;

				X1 = normalizedX1 * (double)WINDOW_WIDTH;
				X2 = normalizedX2 * (double)WINDOW_WIDTH;
				X3 = normalizedX3 * (double)WINDOW_WIDTH;

				Y1 = normalizedY1 * (double)WINDOW_HEIGHT;
				Y2 = normalizedY2 * (double)WINDOW_HEIGHT;
				Y3 = normalizedY2 * (double)WINDOW_HEIGHT;

				HSV1 = normalizedX1 * 360.0;
				HSV2 = normalizedX2 * 360.0;
				HSV3 = normalizedX3 * 360.0;

				X1 *= zoom;
				X2 *= zoom;
				X3 *= zoom;

				Y1 *= zoom;
				Y2 *= zoom;
				Y3 *= zoom;

				X1 -= actShiftX;
				X2 -= actShiftX;
				X3 -= actShiftX;

				Y1 -= actShiftY;
				Y2 -= actShiftY;
				Y3 -= actShiftY;

				//ottimizzazione, disegno ciò che vedo, si può fare anche prima però ci devo pensare!!!
				if (( X1 < 0 || X1 > WINDOW_WIDTH) && (X2 < 0 || X2 > WINDOW_WIDTH) && (X3 < 0 || X3 > WINDOW_WIDTH))
					continue;
				else if((Y1 < 0 || Y1 > WINDOW_HEIGHT) && (Y2 < 0 || Y2 > WINDOW_HEIGHT) && (Y3 < 0 || Y3 > WINDOW_HEIGHT))
					continue;

				v1.position = sf::Vector2f(X1, Y1);
				v2.position = sf::Vector2f(X2, Y2);
				v3.position = sf::Vector2f(X3, Y3);

				if (0)
				{
					sf::Uint8 r1 = HSVTable[(int)HSV1].r * 255;
					sf::Uint8 g1 = HSVTable[(int)HSV1].g * 255;
					sf::Uint8 b1 = HSVTable[(int)HSV1].b * 255;

					sf::Uint8 r2 = HSVTable[(int)HSV2].r * 255;
					sf::Uint8 g2 = HSVTable[(int)HSV2].g * 255;
					sf::Uint8 b2 = HSVTable[(int)HSV2].b * 255;

					sf::Uint8 r3 = HSVTable[(int)HSV3].r * 255;
					sf::Uint8 g3 = HSVTable[(int)HSV3].g * 255;
					sf::Uint8 b3 = HSVTable[(int)HSV3].b * 255;

					v1.color = sf::Color(r1, g1, b1);
					v2.color = sf::Color(r2, g2, b2);
					v3.color = sf::Color(r3, g3, b3);
				}
				v1.color = sf::Color::Red;
				v2.color = sf::Color::Green;
				v3.color = sf::Color::Blue;

				if (0)
				{
					lineA[0].position = sf::Vector2f(X1, Y1);
					lineA[0].color = sf::Color::White;
					lineA[1].position = sf::Vector2f(X2, Y2);
					lineA[1].color = sf::Color::White;

					lineB[0].position = sf::Vector2f(X2, Y2);
					lineB[0].color = sf::Color::White;
					lineB[1].position = sf::Vector2f(X3, Y3);
					lineB[1].color = sf::Color::White;

					lineC[0].position = sf::Vector2f(X3, Y3);
					lineC[0].color = sf::Color::White;
					lineC[1].position = sf::Vector2f(X1, Y1);
					lineC[1].color = sf::Color::White;
				}
#pragma omp critical
				{
					trianglesToDraw.append(v1);
					trianglesToDraw.append(v2);
					trianglesToDraw.append(v3);
				}
			}
			drawTriangles = FALSE;
		}
	}

	return 0;
}

int main()
{
	std::string folder(EXAMPLE_FOLDER);

	printf("nome input: ");
	std::cin >> inputName;

	inputName = folder + inputName;

	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;

	char titleWindow[128];
	double fps = 0;
	bool open = true;

	triangulatioDone = FALSE;
	isRunning = TRUE;

	for (int i = 0; i < 360; i++)
	{
		HSVTable[i] = hsv2rgb({(double)i, 1, 1});
	}

	sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Triangolazione");

	//window.setActive(false);

	ImGui::SFML::Init(window);

	threadTriangulation = CreateThread(
		NULL,                   // default security attributes
		0,                      // use default stack size  
		computeTriangulation,       // thread function name
		NULL,          // argument to thread function 
		0,                      // use default creation flags 
		&threadTrID);   // returns the thread identifier 

	threadActions = CreateThread(
		NULL,                   // default security attributes
		0,                      // use default stack size  
		actionThread,       // thread function name
		NULL,          // argument to thread function 
		0,                      // use default creation flags 
		&threadActID);   // returns the thread identifier 

	//window.setActive(true);

	bool single = false;
	sf::Clock deltaClock;
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			ImGui::SFML::ProcessEvent(event);
			// "close requested" event: we close the window
			if (event.type == sf::Event::Closed)
				window.close();

			//if (event.type == sf::Event::MouseButtonPressed)
			//{
			//	window.clear();

			//	single = !single;

			//	if (single)
			//	{
			//		DrawTriangles(&window, xMin, xMax, yMin, yMax);

			//		window.display();
			//	}
			//	else
			//	{
			//		DrawTrianglesThreads(&window, xMin, xMax, yMin, yMax);

			//		window.display();
			//	}
			//}
		}

		ImGui::SFML::Update(window, deltaClock.restart());

		sprintf_s(titleWindow, "fps: %f", fps);

		sf::String title(titleWindow);

		start = std::chrono::high_resolution_clock::now();
		
		window.setTitle(title);

		ImGui::Begin("parametri triangoli", &open);
		ImGui::SliderFloat("Zoom", &zoom, 1.0, 200.0);
		ImGui::SliderFloat("Shift X", &shiftX, 0, WINDOW_WIDTH);
		ImGui::SliderFloat("Shift Y", &shiftY, 0, WINDOW_HEIGHT);
		if(ImGui::Button("Applica Impostazioni"))
		{
			if(triangulatioDone)
				drawTriangles = TRUE;
		}
		ImGui::End();

		window.clear(sf::Color::Black);
		if (triangulatioDone && !drawTriangles)
		{
			window.draw(trianglesToDraw);
			//DrawTriangles(&window, xMin, xMax, yMin, yMax);
			//DrawTrianglesThreads(&window, xMin, xMax, yMin, yMax);
		}

		//ImGui::ShowDemoWindow();
		ImGui::SFML::Render(window);

		window.display();

		end = std::chrono::high_resolution_clock::now();

		fps = (float)1e9 / (float)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	}

	for (int i = 0; i < NUM_THREADS; i++)
	{
		delete triangulators[i];
	}

	delete triangulator;

	ImGui::SFML::Shutdown();
	isRunning = FALSE;

	return 0;
}