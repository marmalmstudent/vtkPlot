#include "Plot.hpp"

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <array>

#include "root_directory.h"

union cuf {
  char c[sizeof(float)];
  float f;
};

std::vector<float>
Plot::BytesToDouble(std::vector<char> const &src)
{
  std::vector<float> dst;
  std::vector<char>::const_iterator itr = src.begin();
  for (; itr != src.end();) {
    cuf uni;
    char *c = uni.c;
    for(size_t j = sizeof(float); j > 0 && itr != src.end(); --j) {
      *c++ = *itr++;
    }
    dst.push_back(uni.f);
  }
  return dst;
}

std::vector<char>
Plot::ReadAllBytes(std::string const &filename)
{
  std::vector<char> result;
  std::ifstream ifs(filename, std::ifstream::binary | std::ifstream::ate);
  if (ifs.is_open()) {
    std::ifstream::pos_type pos = ifs.tellg();
    result.resize(pos);
    ifs.seekg(0, std::ios::beg);
    ifs.read(&result[0], pos);
    ifs.close();
  } else
    cout << "unable to find file: " << filename << endl;
  return result;
}

Bounds<float>
Plot::calcZBounds(std::vector<float> const &arr)
{
  Bounds<float> zBounds = { arr[0], arr[0] };
  for (float const &f : arr) {
    if (f > zBounds.max) {
      zBounds.max = f;
    }
    if (f < zBounds.min) {
      zBounds.min = f;
    }
  }
  return zBounds;
}

Bounds<float>
Plot::createTopology(vtkSmartPointer<vtkPolyData> point,
		     std::vector<float> const &x,
		     std::vector<float> const &y,
		     std::vector<float> const &z)
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

vtkSmartPointer<vtkPolyDataMapper>
Plot::createSurface()
{
  vtkSmartPointer<vtkPolyDataMapper> surface = vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkPolyData> points = vtkSmartPointer<vtkPolyData>::New();
  Bounds<float> zBounds = createTopology(points, loadData(fileX), loadData(fileY), loadData(fileZ));
  surface->SetInputData(points);
  surface->SetScalarRange(zBounds.min, zBounds.max); // colormap boundaries
  return surface;
}

vtkSmartPointer<vtkActor>
Plot::makeActor(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetOpacity(1.0);
  actor->RotateX(-45);
  actor->RotateY(45);
  return actor;
}

vtkSmartPointer<vtkRenderer>
Plot::makeRender(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
  vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
  ren->AddActor(makeActor(mapper)); // add the rendering scene
  ren->SetBackground(0.7, 0.8, 1.0);
  return ren;
}

vtkSmartPointer<vtkRenderWindow>
Plot::makeRenderWindow(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(makeRender(mapper)); // add renderer to window
  renWin->SetSize(800, 600);
  return renWin;
}

void
Plot::plotSurface(vtkSmartPointer<vtkPolyDataMapper> mapper)
{
  // Interactor to interact with the render window
  vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(makeRenderWindow(mapper));
  iren->Render();
  iren->Start();
}

int
main(int argc,
     char *argv[])
{
  if (argc > 5) {
    std::string fileX(sys_path(argv[1]));
    std::string fileY(sys_path(argv[2]));
    std::string fileZ(sys_path(argv[3]));
    int nxpts = std::stoll(argv[4]);
    int nypts = std::stoll(argv[5]);
    
    std::cout << "Reading from file:";
    std::cout << "\n\t" << fileX;
    std::cout << "\n\t" << fileY;
    std::cout << "\n\t" << fileZ << std::endl;
    
    std::cout << "Grid size: " << nxpts << " x " << nypts << std::endl;
    
    Plot p(fileX, fileY, fileZ, nxpts, nypts);
    p.plotSurface(p.createSurface());
    
    return EXIT_SUCCESS;
  } else {
    std::cout << "Usage: " << std::string(argv[0]) << " <xfile> <yfile> <zfile> <width> <height>" << std::endl;
    return EXIT_FAILURE;
  }
}
