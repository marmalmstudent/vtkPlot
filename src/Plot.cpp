#include "Plot.hpp"

#include <vtk/vtkVersion.h>
#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkFloatArray.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkActor.h>
#include <vtk/vtkProperty.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderer.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkCamera.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <array>

using namespace std;

union cuf {
	char c[sizeof(float)];
	float f;
};

vector<float> Plot::BytesToDouble(const vector<char> &src)
{
	vector<float> dst;
	for (unsigned i = 0; i < src.size(); i += sizeof(float)) {
		cuf uni;
		for(unsigned j = 0; j < sizeof(float); j++)
			uni.c[j] = src[i+j];
		dst.push_back(uni.f);
	}
	return dst;
}

vector<char> Plot::ReadAllBytes(const string &filename)
{
	vector<char> result;
	ifstream ifs(filename, ifstream::binary | ifstream::ate);
	if (ifs.is_open()) {
		ifstream::pos_type pos = ifs.tellg();
		result.resize(pos);
		ifs.seekg(0, ios::beg);
		ifs.read(&result[0], pos);
		ifs.close();
	} else
		cout << "unable to find file: " << filename << endl;
	return result;
}

array<float, 2> Plot::calcZBounds(const vector<float> &arr)
{
	array<float, 2> zBounds;
	zBounds[0] = arr[0];
	zBounds[1] = arr[0];
	for (vector<float>::const_iterator itr = arr.begin(); itr != arr.end(); ++itr) {
		if (*itr > zBounds[1])
			zBounds[1] = *itr;
		if (*itr < zBounds[0])
			zBounds[0] = *itr;
    }
	return zBounds;
}

array<float, 2> Plot::createTopology(vtkSmartPointer<vtkPolyData> point,
		const vector<float> &x, const vector<float> &y,
		const vector<float> &z)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> colormap = vtkSmartPointer<vtkFloatArray>::New();

	for (int i = 0; i < xSize-1; ++i) {
		for (int j = 0; j < ySize-1; ++j) {
			vtkIdType p[] = {
					points->InsertNextPoint(x[(j  )*xSize + (i  )], y[(j  )*xSize + (i  )], z[(j  )*xSize + (i  )]),
					points->InsertNextPoint(x[(j  )*xSize + (i+1)], y[(j  )*xSize + (i+1)], z[(j  )*xSize + (i+1)]),
					points->InsertNextPoint(x[(j+1)*xSize + (i  )], y[(j+1)*xSize + (i  )], z[(j+1)*xSize + (i  )]),
					points->InsertNextPoint(x[(j+1)*xSize + (i+1)], y[(j+1)*xSize + (i+1)], z[(j+1)*xSize + (i+1)])
			};
			vertices->InsertNextCell(4, p);
			colormap->InsertNextValue(z[(j  )*xSize + (i  )]);
			colormap->InsertNextValue(z[(j  )*xSize + (i+1)]);
			colormap->InsertNextValue(z[(j+1)*xSize + (i  )]);
			colormap->InsertNextValue(z[(j+1)*xSize + (i+1)]);
		}
    }
	point->SetPoints(points);
	//point->SetVerts(vertices); // points, slower than surface
	point->SetStrips(vertices); // surface, faster than points
	point->GetPointData()->SetScalars(colormap);
	return calcZBounds(z);
}

vtkSmartPointer<vtkPolyDataMapper> Plot::createSurface()
{
	vtkSmartPointer<vtkPolyDataMapper> surface = vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkPolyData> points = vtkSmartPointer<vtkPolyData>::New();
	array<float, 2> zBounds = createTopology(points, loadData(fileX), loadData(fileY), loadData(fileZ));
	surface->SetInputData(points);
	surface->SetScalarRange(zBounds[0], zBounds[1]); // colormap boundaries
	return surface;
}

vtkSmartPointer<vtkActor> Plot::makeActor(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetOpacity(1.0);
	actor->RotateX(-45);
	actor->RotateY(45);
	return actor;
}

vtkSmartPointer<vtkRenderer> Plot::makeRender(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	ren->AddActor(makeActor(mapper)); // add the rendering scene
	ren->SetBackground(0.7, 0.8, 1.0);
	return ren;
}

vtkSmartPointer<vtkRenderWindow> Plot::makeRenderWindow(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(makeRender(mapper)); // add renderer to window
	renWin->SetSize(800, 600);
	return renWin;
}

void Plot::plotSurface(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
	// Interactor to interact with the render window
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(makeRenderWindow(mapper));
	iren->Render();
	iren->Start();
}

int main(int argc, char **argv)
{
	if (argc > 5) {
		string fileX = argv[1];
		string fileY = argv[2];
		string fileZ = argv[3];
		int nxpts = stoi(argv[4]);
		int nypts = stoi(argv[5]);
		Plot p(fileX, fileY, fileZ, nxpts, nypts);
		printf("Reading from file: \n\t%s\n\t%s\n\t%s\n", fileX.c_str(), fileY.c_str(), fileZ.c_str());
		printf("Grid size: %d x %d", nxpts, nypts);
		p.plotSurface(p.createSurface());
		return EXIT_SUCCESS;
	} else {
		printf("Usage: %s <xfile> <yfile> <zfile> <width> <height>", argv[0]);
		return EXIT_FAILURE;
	}
}
