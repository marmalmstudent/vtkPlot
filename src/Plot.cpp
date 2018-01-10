/*
  This code is optimized for 80 columns
|------------------------------------------------------------------------------|
 */
#include <vtk/vtkVersion.h>
#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkDoubleArray.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkActor.h>
#include <vtk/vtkProperty.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderer.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkCamera.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>  //for std::istringstream
#include <iterator> //for std::istream_iterator
#include <vector>   //for std::vector

#define VTK_SP(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
#define PI 3.141592653589793238463

using namespace std;

union cuf {
	char c[sizeof(float)];
	float f;
};

void BytesToDouble(vector<char> &src, vector<float> &dst, bool littleEndian)
{
	for (unsigned i = 0; i < src.size(); i += sizeof(float)) {
		cuf uni;
		for(unsigned j = 0; j < sizeof(float); j++)
			uni.c[j] = src[i+j];
		dst.push_back(uni.f);
	}
}

static void ReadAllBytes(string &filename, vector<char> &result)
{
	ifstream ifs(filename, ifstream::binary | ifstream::ate);
	if (ifs.is_open()) {
		ifstream::pos_type pos = ifs.tellg();
		result.resize(pos);

		ifs.seekg(0, ios::beg);
		ifs.read(&result[0], pos);
		ifs.close();
	} else {
		cout << "unable to find file: " << filename << endl;
	}
}

void loadData(string &filename, vector<float> &arr)
{
  vector<char> b;
  ReadAllBytes(filename, b);
  BytesToDouble(b, arr, true);
}

void calcZBounds(vector<float> &arr, vector<float> *zBounds)
{
	for (vector<float>::iterator itr = arr.begin(); itr != arr.end(); ++itr) {
		if (*itr > zBounds->at(1))
			zBounds->at(1) = *itr;
		if (*itr < zBounds->at(0))
			zBounds->at(0) = *itr;
    }
}

void createTopology(vtkSmartPointer<vtkPoints> &points,
                    vtkSmartPointer<vtkCellArray> &vertices,
                    vtkSmartPointer<vtkDoubleArray> &colormap,
                    vector<float> &x, vector<float> &y, vector<float> &z,
                    int xSize, int ySize)
{
	vtkIdType pid[4];
	for (int i = 0; i < xSize-1; ++i) {
		for (int j = 0; j < ySize-1; ++j) {
			// The coordinates of the points in the cell
			// Insert the points
			pid[0] = points->InsertNextPoint(x[j*xSize+i], y[j*xSize+i], z[j*xSize+i]);
			pid[1] = points->InsertNextPoint(x[j*xSize+i+1], y[j*xSize+i+1], z[j*xSize+i+1]);
			pid[2] = points->InsertNextPoint(x[(j+1)*xSize+i], y[(j+1)*xSize+i], z[(j+1)*xSize+i]);
			pid[3] = points->InsertNextPoint(x[(j+1)*xSize+i+1], y[(j+1)*xSize+i+1], z[(j+1)*xSize+i+1]);
			// Insert the next cell
			vertices->InsertNextCell(4, pid);
			// Insert height to colormap
			colormap->InsertNextValue(z[j*xSize +i]);
			colormap->InsertNextValue(z[j*xSize +i+1]);
			colormap->InsertNextValue(z[(j+1)*xSize +i]);
			colormap->InsertNextValue(z[(j+1)*xSize +i+1]);
		}
    }
}

void plotSurface(vtkRenderWindowInteractor &iren, string &fileX, string &fileY, string &fileZ, int xSize, int ySize)
{
	/* Create the geometry of a point (the coordinate) */
	VTK_SP(vtkPoints, points);
	/* Create the topology of the point (a vertex) */
	VTK_SP(vtkCellArray, vertices);
	/* Create the color of the point based on topology */
	VTK_SP(vtkDoubleArray, colormap);

	/*
    ----------------------------------------------------------------------------
    z-coords need unique values for every set of x and y. x and y are constant
    along y and x respecively so no they don't have to be same size as z.
    Allocate in the heap to manage memory.
    ----------------------------------------------------------------------------
	 */
	/* Create the surface */
	vector<float> x, y, z;
	loadData(fileX, x);
	loadData(fileY, y);
	loadData(fileZ, z);
	/* Find max and min for the colormap */
	vector<float> zBounds;
	zBounds.push_back(z[0]);
	zBounds.push_back(z[0]);
	calcZBounds(z, &zBounds);

	/* Create the topology of the point (a vertex) */
	createTopology(points, vertices, colormap, x, y, z, xSize, ySize);
	/* Clean up unused variables */
	x.clear(); x.shrink_to_fit();
	y.clear(); y.shrink_to_fit();
	z.clear(); z.shrink_to_fit();

	/* Create a polydata object */
	VTK_SP(vtkPolyData, point);
	/* Map polygon data to graphics */
	VTK_SP(vtkPolyDataMapper, mapper);
	/* Represent graphics in a rendering scene */
	VTK_SP(vtkActor, actor);

	/* create polygons data from then points */
	point->SetPoints(points);
	//point->SetVerts(vertices); // points, slower than surface
	point->SetStrips(vertices); // surface, faster than points
	point->GetPointData()->SetScalars(colormap);

	/* Visualize polygons as primitive grapnics */
	mapper->SetInputData(point);
	mapper->SetScalarRange(zBounds[0], zBounds[1]); // colormap boundaries

	/* set the rendering scene */
	actor->SetMapper(mapper);
	actor->GetProperty()->SetOpacity(1.0);
	actor->RotateX(-45);
	actor->RotateY(45);

	/*
    ----------------------------------------------------------------------------
    Add then rendering scene to the renderer and add it to the rendering window.
    Finally add an interactor to enable interaction with the rendering window.
    ----------------------------------------------------------------------------
	 */
	/* Controls rendering of objects */
	VTK_SP(vtkRenderer, ren);
	/* Window to place the renderer in */
	VTK_SP(vtkRenderWindow, renWin);

	ren->AddActor(actor); // add the rendering scene
	ren->SetBackground(0.7, 0.8, 1.0);
	renWin->SetSize(800, 800);
	renWin->AddRenderer(ren); // add renderer to window
	iren.SetRenderWindow(renWin); // set the window for the interactor
}

int main(int argc, char **argv)
{
	int nxpts, nypts;
	string fileX, fileY, fileZ;
	if (argc > 5) {
		fileX = argv[1];
		fileY = argv[2];
		fileZ = argv[3];
		nxpts = stoi(argv[4]);
		nypts = stoi(argv[5]);
	} else {
		printf("Usage: %s <xfile> <yfile> <zfile> <width> <height>", argv[0]);
		return EXIT_FAILURE;
	}
	printf("Reading from file: \n\t%s\n\t%s\n\t%s\n", fileX.c_str(), fileY.c_str(), fileZ.c_str());
	printf("Grid size: %d x %d", nxpts, nypts);
	// Interactor to interact with the render window
	VTK_SP(vtkRenderWindowInteractor, iren);
	//std::vector<double> surfData = loadData(filename, nxpts, nypts);
	plotSurface(*iren, fileX, fileY, fileZ, nxpts, nypts);
	//renWin->Render();
	//ren->GetActiveCamera()->Zoom(1.5);
	iren->Render();
	iren->Start();
	return EXIT_SUCCESS;
}
