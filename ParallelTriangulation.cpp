// ParallelTriangulation.cpp : Questo file contiene la funzione 'main', in cui inizia e termina l'esecuzione del programma.
//
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <CDT.h>
#include <omp.h>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <Windows.h>
#include <chrono>
#include <unordered_map>

#include "imgui.h"
#include "imgui-SFML.h"

//strutture dati custom
#include "mainutility.h"

#pragma comment(lib, "sfml-system.lib")
#pragma comment(lib, "sfml-graphics.lib")
#pragma comment(lib, "sfml-window.lib")
#pragma comment(lib, "opengl32.lib")

#define NUM_THREADS 12
#define EXAMPLE_FOLDER "examples/"
#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600

CDT::Triangulation<double> triangulator; //singolo thread
std::vector<CDT::V2d<double>> points;
std::vector<double> zHeight;

std::vector<double> zHeights[NUM_THREADS];
std::vector<CDT::V2d<double>> tPoints[NUM_THREADS];
CDT::Triangulation<double>* triangulators[NUM_THREADS];
Quadrante quadranti[NUM_THREADS];

std::unordered_map<CDT::V2d<double>, bool> hashTablePoints; //per togliere duplicati vari

std::vector<XYZ> surface;

std::string inputName;

sf::VertexArray trianglesToDraw(sf::Triangles);
sf::VertexArray trianglesPToDraw(sf::Triangles);

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
bool drawTrianglesP; //con computo parallelo
bool computeSurface;
bool saveOBJ;
bool isRunning;

float zoom = 1.0;
float shiftX = 0;
float shiftY = 0;
float gridResolution = 0;

rgb HSVTable[360];
std::string folder;

double start, end;

void createOBJ()
{
	FILE* fp;
	fopen_s(&fp, "mesh.obj", "wb");

	if (fp)
	{
		for (int i = 0; i < triangulator.vertices.size(); i++)
		{
			fprintf(fp, "v %f %f %f\n", triangulator.vertices[i].x, triangulator.vertices[i].y, zHeight[i]);
		}

		for (int i = 0; i < triangulator.triangles.size(); i++)
		{
			fprintf(fp, "f %d %d %d\n", 
				triangulator.triangles[i].vertices[0] + 1, 
				triangulator.triangles[i].vertices[1] + 1,
				triangulator.triangles[i].vertices[2] + 1);
		}

		fclose(fp);
	}

	fopen_s(&fp, "mesh_parallel.obj", "wb");

	if (fp)
	{
		std::vector<int> offset; //sennò gli indici non sono corretti

		for (int j = 0; j < NUM_THREADS; j++)
		{
			for (int i = 0; i < triangulators[j]->vertices.size(); i++)
			{
				fprintf(fp, "v %f %f %f\n", triangulators[j]->vertices[i].x, triangulators[j]->vertices[i].y, zHeights[j][i]);
			}

			if (j == 0)
				offset.push_back(0);
			else
				offset.push_back(offset[j - 1] + triangulators[j - 1]->vertices.size());
		}

		for (int j = 0; j < NUM_THREADS; j++)
		{
			for (int i = 0; i < triangulators[j]->triangles.size(); i++)
			{
				fprintf(fp, "f %d %d %d\n",
					triangulators[j]->triangles[i].vertices[0] + 1 + offset[j],
					triangulators[j]->triangles[i].vertices[1] + 1 + offset[j],
					triangulators[j]->triangles[i].vertices[2] + 1 + offset[j]);
			}
		}

		fclose(fp);
	}
}

void saveSurface()
{
	FILE* fp;
	fopen_s(&fp, "surface.xyz", "wb");

	if (fp)
	{
		for (int i = 0; i < surface.size(); i++)
		{
			fprintf(fp, "%f %f %f\n", surface[i].x, surface[i].y, surface[i].z);
		}

		fclose(fp);
	}
}

void saveSurface(std::vector<CDT::V2d<double>>& XY, std::vector<double>& Z, std::vector<int>& mask)
{
	FILE* fp;
	fopen_s(&fp, "surface.xyz", "wb");

	if (fp)
	{
		for (int i = 0; i < XY.size(); i++)
		{
			if (mask[i] > 0)
			{
				fprintf(fp, "%f %f %f\n", XY[i].x, XY[i].y, Z[i]);
			}
		}

		fclose(fp);
	}
}

DWORD WINAPI computeTriangulation(LPVOID lpParam)
{
	int npoints = 0;
	double x, y, z;

	std::cout << "sizeof(triangulator) in byte: " << sizeof(triangulator) << std::endl;

	FILE* fp;
	
	{
		do
		{
			printf("nome input: ");
			std::cin >> inputName;

			inputName = folder + inputName;

			fopen_s(&fp, inputName.c_str(), "rb");
		}
		while (fp == NULL);

		printf("leggo da file\n");

		fscanf_s(fp, "%d\n", &npoints); //prima riga deve contenere numero punti

		printf("numero punti: %d\n", npoints);

		CDT::V2d<double> newPoint;

		for (int i = 0; i < npoints; i++)
		{
			fscanf_s(fp, "%lf %lf %lf \n", &x, &y, &z);

			newPoint = { x, y };

			if (hashTablePoints[newPoint])
				continue;
			else
				hashTablePoints[newPoint] = true;

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

		npoints = points.size(); //riprendo la dimensione vera dopo il filtraggio

		printf("numero punti dopo filtraggio: %d\n", npoints);
	}

	{
		printf("separo i punti e rimuovo duplicati\n");

		double latoQuadrante = abs(xMax - xMin) / (double)NUM_THREADS;

		for (int i = 0; i < NUM_THREADS; i++)
		{
			double start = xMin + latoQuadrante * i;
			double end = start + latoQuadrante; //aggiungo un pò di sovrapposizione

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

		//triangolazioni
		if(1)
		{
			printf("inizio single thread\n");
			start = omp_get_wtime();
			{
				triangulator.insertVertices(points);
				triangulator.eraseSuperTriangle();
			}
			end = omp_get_wtime();

			std::cout << "sequenziale tempo: " << end - start << " sec" << std::endl;

			printf("inizio parallelo\n");

			start = omp_get_wtime();

#pragma omp parallel num_threads(NUM_THREADS) shared(triangulators, tPoints)
			{
				int curThread = omp_get_thread_num();

				triangulators[curThread]->insertVertices(tPoints[curThread]);
				triangulators[curThread]->eraseSuperTriangle();
			}

			end = omp_get_wtime();

			std::cout << "parallelo: " << end - start << " sec" << std::endl;
		}
	}

	drawTriangles = TRUE;
	while (drawTriangles)
	{
		Sleep(1);
	}

	drawTrianglesP = TRUE;
	while (drawTrianglesP)
	{
		Sleep(1);
	}

	triangulatioDone = TRUE;

	return 0;
}

DWORD WINAPI actionThread(LPVOID lpParam)
{
	while (isRunning)
	{
		if (drawTriangles)
		{
			//if (triangulator == NULL)
			/*{
				drawTriangles = FALSE;
				continue;
			}*/

			std::vector<CDT::V2d<double>>& curPoints = points;

			auto start = omp_get_wtime();

			double xDiff = 1.0 / (xMax - xMin);
			double yDiff = 1.0 / (yMax - yMin);
			trianglesToDraw.clear();
			double actShiftX = shiftX * zoom;
			double actShiftY = shiftY * zoom;

			int numTriangles = triangulator.triangles.size();

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

				CDT::Triangle* triangle = &triangulator.triangles[j];

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

				if (1)
				{
					sf::Color tintaV1, tintaV2, tintaV3;
					tintaV1.r = (BYTE)(HSVTable[(int)HSV1].r * 255);
					tintaV1.g = (BYTE)(HSVTable[(int)HSV1].g * 255);
					tintaV1.b = (BYTE)(HSVTable[(int)HSV1].b * 255);

					tintaV2.r = (BYTE)(HSVTable[(int)HSV2].r * 255);
					tintaV2.g = (BYTE)(HSVTable[(int)HSV2].g * 255);
					tintaV2.b = (BYTE)(HSVTable[(int)HSV2].b * 255);

					tintaV3.r = (BYTE)(HSVTable[(int)HSV3].r * 255);
					tintaV3.g = (BYTE)(HSVTable[(int)HSV3].g * 255);
					tintaV3.b = (BYTE)(HSVTable[(int)HSV3].b * 255);

					v1.color = tintaV1;
					v2.color = tintaV2;
					v3.color = tintaV3;
				}

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
#pragma omp critical //da cambiare
				{
					trianglesToDraw.append(v1);
					trianglesToDraw.append(v2);
					trianglesToDraw.append(v3);
				}
			}

			auto end = omp_get_wtime();

			std::cout << "tempo disegno triangulator singolo: " << end - start << " sec" << std::endl;

			drawTriangles = FALSE;
		}
		if (drawTrianglesP)
		{
			auto start = omp_get_wtime();

			double xDiff = 1.0 / (xMax - xMin);
			double yDiff = 1.0 / (yMax - yMin);
			trianglesPToDraw.clear();
			double actShiftX = shiftX * zoom;
			double actShiftY = shiftY * zoom;
			
#pragma omp parallel num_threads(NUM_THREADS)
			{
				int curTr = omp_get_thread_num();

				int numTriangles = triangulators[curTr]->triangles.size();

				std::vector<CDT::V2d<double>>& curPoints = tPoints[curTr];

				for (int j = 0; j < numTriangles; j++)
				{
					sf::Vertex v1;
					sf::Vertex v2;
					sf::Vertex v3;

					CDT::Triangle* triangle = &triangulators[curTr]->triangles[j];

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
					if ((X1 < 0 || X1 > WINDOW_WIDTH) && (X2 < 0 || X2 > WINDOW_WIDTH) && (X3 < 0 || X3 > WINDOW_WIDTH))
						continue;
					else if ((Y1 < 0 || Y1 > WINDOW_HEIGHT) && (Y2 < 0 || Y2 > WINDOW_HEIGHT) && (Y3 < 0 || Y3 > WINDOW_HEIGHT))
						continue;

					v1.position = sf::Vector2f(X1, Y1);
					v2.position = sf::Vector2f(X2, Y2);
					v3.position = sf::Vector2f(X3, Y3);

					if (1)
					{
						sf::Color tintaV1, tintaV2, tintaV3;
						tintaV1.r = (BYTE)(HSVTable[(int)HSV1].r * 255);
						tintaV1.g = (BYTE)(HSVTable[(int)HSV1].g * 255);
						tintaV1.b = (BYTE)(HSVTable[(int)HSV1].b * 255);

						tintaV2.r = (BYTE)(HSVTable[(int)HSV2].r * 255);
						tintaV2.g = (BYTE)(HSVTable[(int)HSV2].g * 255);
						tintaV2.b = (BYTE)(HSVTable[(int)HSV2].b * 255);

						tintaV3.r = (BYTE)(HSVTable[(int)HSV3].r * 255);
						tintaV3.g = (BYTE)(HSVTable[(int)HSV3].g * 255);
						tintaV3.b = (BYTE)(HSVTable[(int)HSV3].b * 255);

						v1.color = tintaV1;
						v2.color = tintaV2;
						v3.color = tintaV3;
					}

#pragma omp critical //da cambiare
					{
						trianglesPToDraw.append(v1);
						trianglesPToDraw.append(v2);
						trianglesPToDraw.append(v3);
					}
				}
			}
			auto end = omp_get_wtime();

			std::cout << "tempo disegno triangulators threads: " << end - start << " sec" << std::endl;

			drawTrianglesP = FALSE;
		}
		if (computeSurface)
		{
			std::cout << "inizio computo griglia parallelo" << std::endl;

			//metodo artigianale fatto da simone...
			int numX = 0;
			int numY = 0;

			for (double j = yMin; j <= yMax; j += gridResolution)
				numY++;

			for (double i = xMin; i <= xMax; i += gridResolution)
				numX++;

			start = omp_get_wtime();
			std::vector<CDT::V2d<double>> xySteps;
			std::vector<double> surfHeight; //altezza superficie
			std::vector<int> surfHeightCount; //numero campioni superficie

			//com'è fatta la griglia
			/*
				(X1, Y1), (X2, Y1) .....
				(X1, Y2), (X2, Y2) .....
				.....	  .....
			*/

			for (double j = yMin; j <= yMax; j += gridResolution)
			{
				for (double i = xMin; i <= xMax; i += gridResolution)
				{
					xySteps.push_back({ i, j });
					surfHeight.push_back(0); //piatta inizialmente (pianerottolo...)
					surfHeightCount.push_back(0);
				}
			}

			double stepGrid = gridResolution / 2.0;

#pragma omp parallel for
			for (int i = 0; i < points.size(); i++)
			{
				int tr = omp_get_thread_num();

				//double bestX = abs(points[i].x - xySteps[0].x);
				//int indX = 0;

				//double bestY = abs(points[i].y - xySteps[0].y);
				//int indY = 0;

				//yCandidates[tr].arr.clear();
				//xCandidates[tr].arr.clear();

				int xCandidate = 0;
				int yCandidate = 0;

				double xminP = points[i].x - stepGrid;
				double xmaxP = points[i].x + stepGrid;
				double yminP = points[i].y - stepGrid;
				double ymaxP = points[i].y + stepGrid;

				for (int j = 0; j < numX; j++)
				{
					if (xminP <= xySteps[j].x && xySteps[j].x <= xmaxP)
					{
						xCandidate = j;
						break;
						//xCandidates[tr].arr.push_back(j);
					}
				}

				for (int j = 0; j < numY; j++)
				{
					if (yminP <= xySteps[j * numX].y && xySteps[j * numX].y <= ymaxP)
					{
						yCandidate = j * numX;
						break;
						//yCandidates[tr].arr.push_back(j);
					}
				}

#pragma omp critical
				{
					surfHeight[yCandidate + xCandidate] += zHeight[i];
					surfHeightCount[yCandidate + xCandidate]++;
				}
			}

#pragma omp parallel for
			for (int i = 0; i < surfHeight.size(); i++)
			{
				if (surfHeightCount[i])
					surfHeight[i] = surfHeight[i] / (double)surfHeightCount[i];
			}

			end = omp_get_wtime();
			std::cout << "metodo artigianale griglia parallelo: " << end - start << " sec" << std::endl;

			std::cout << "salvataggio su file (surface.xyz)" << std::endl;
			saveSurface(xySteps, surfHeight, surfHeightCount);
			
			computeSurface = FALSE;
		}
		if (saveOBJ)
		{
			std::cout << "inizio obj mesh (salvataggio su file)" << std::endl;
			start = omp_get_wtime();

			createOBJ();

			end = omp_get_wtime();

			std::cout << "creazione obj mesh (salvataggio su file): " << end - start << " sec" << std::endl;

			saveOBJ = FALSE;
		}
	}

	return 0;
}

int main()
{
	computeSurface = FALSE;

	for (int i = 0; i < NUM_THREADS; i++)
	{
		triangulators[i] = 0;
		triangulators[i] = new CDT::Triangulation<double>();
	}

	folder = EXAMPLE_FOLDER;

	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;

	char titleWindow[128];
	double fps = 0;
	bool open = true;
	bool parallelDraw = true;

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
		}

		ImGui::SFML::Update(window, deltaClock.restart());

		sprintf_s(titleWindow, "fps: %f", fps);

		sf::String title(titleWindow);

		start = std::chrono::high_resolution_clock::now();
		
		window.setTitle(title);

		ImGui::Begin("parametri triangoli", &open);
		if (!triangulatioDone)
		{
			ImGui::Text("Computando l'input nuvola di punti........");
		}
		ImGui::SliderFloat("Zoom", &zoom, 1.0, 200.0);
		ImGui::SliderFloat("Shift X", &shiftX, 0, WINDOW_WIDTH);
		ImGui::SliderFloat("Shift Y", &shiftY, 0, WINDOW_HEIGHT);
		if(ImGui::Button("Applica Impostazioni"))
		{
			if (triangulatioDone)
				drawTriangles = TRUE;

			if (triangulatioDone)
				drawTrianglesP = TRUE;
		}
		ImGui::Checkbox("Parallelo", &parallelDraw);
		ImGui::SliderFloat("Risoluzione griglia superfice", &gridResolution, 0.1, 2);
		if (ImGui::Button("Computa Superfice in griglia"))
		{
			if(!computeSurface)
				computeSurface = TRUE;
		}
		if (ImGui::Button("Salva mesh"))
		{
			if (triangulatioDone)
				if (!saveOBJ)
					saveOBJ = TRUE;
		}
		ImGui::End();

		window.clear(sf::Color::Black);

		if (parallelDraw)
		{
			if (triangulatioDone && !drawTrianglesP)
			{
				window.draw(trianglesPToDraw);
			}
		}
		else
		{
			if (triangulatioDone && !drawTriangles)
			{
				window.draw(trianglesToDraw);
			}
		}

		ImGui::SFML::Render(window);

		window.display();

		end = std::chrono::high_resolution_clock::now();

		fps = (float)1e9 / (float)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	}

	for (int i = 0; i < NUM_THREADS; i++)
	{
		delete triangulators[i];
	}

	ImGui::SFML::Shutdown();
	isRunning = FALSE;

	return 0;
}