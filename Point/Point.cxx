#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>  //for std::istringstream
#include <iterator> //for std::istream_iterator
#include <vector>   //for std::vector


#define VTK_SP(type, name)\
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#define PI 3.141592653589793238463;

void loadData(std::string filename,
              std::vector<double> *arr)
{
  for (std::vector<double>::iterator i = arr->begin(); i != arr->end(); i++)
    {
      std::cout << *i << std::endl;
    }
  std::string line;
  std::ifstream streamin(filename);
  if (streamin.is_open())
    {
      while(std::getline(streamin, line))
        {
          std::istringstream ss(line);
          std::istream_iterator<std::string> begin(ss), end;
          std::vector<std::string> arrayTokens(begin, end);
          for (int i = 0; i < arrayTokens.size(); i++)
            {
              arr->push_back(std::stod(arrayTokens[i]));
            }
        }
      streamin.close();
    }
  else
    {
      std::cout << "unable to fine file: " << filename << std::endl;
    }
  return;
}

void createSurfaceArrayFromFile(std::vector<double> *x,
                                std::vector<double> *y,
                                std::vector<double> *z,
                                int xSize, int ySize,
                                float xMin, float xMax,
                                float yMin, float yMax,
                                std::string filename)
{
  loadData("/home/marcus/git/vtkPlot/Point/X.txt", x);
  loadData("/home/marcus/git/vtkPlot/Point/Y.txt", y);
  loadData("/home/marcus/git/vtkPlot/Point/Z.txt", z);
}

void createSurfaceArray(float* x, float* y, float* z,
			int xSize, int ySize,
			float xMin, float xMax,
			float yMin, float yMax)
{
  /*
    ----------------------------------------------------------------------------
    Create a surface of sin(x)*cos(y) using evenly spaced values of x and y.
    ----------------------------------------------------------------------------
   */
  float xIncrement = (xMax-xMin)/(xSize-1);
  float yIncrement = (yMax-yMin)/(ySize-1);
  for (int i = 0; i < xSize; ++i)
    {
      for (int j = 0; j < ySize; ++j)
	{
	  x[i] = xMin + i*xIncrement;
	  y[j] = yMin + j*yIncrement;
	  z[j*xSize +i] = sin(x[i])*cos(y[j]);
	}
    }
}

void calcZBounds(std::vector<double> *arr, std::vector<double> *zBounds)
{
  /*
    ----------------------------------------------------------------------------
    Find then maximum and minimum values of a given array to calibrate colormap.
    ----------------------------------------------------------------------------
   */

  for (std::vector<double>::iterator i = arr->begin(); i != arr->end(); i++)
    {
      if (*i > zBounds->at(1))
        zBounds->at(1) = *i;
      if (*i < zBounds->at(0))
        zBounds->at(0) = *i;
    }
}

void createTopology(vtkSmartPointer<vtkPoints> points,
                    vtkSmartPointer<vtkCellArray> vertices,
                    vtkSmartPointer<vtkDoubleArray> colormap,
                    std::vector<double> *x,
                    std::vector<double> *y,
                    std::vector<double> *z,
                    int xSize, int ySize)
{
  /*
    ----------------------------------------------------------------------------
    Add points, create a cell using four points as boundaries and set the color
    based on then elevation (z-value) of the point.
    x,y ---- x+1,y
     |        |
     |  Cell  |
     |        |
    x,y+1 -- x+1,y+1
    ----------------------------------------------------------------------------
   */
  // 4 points per cell
  /*
  vtkIdType pid[4];
  float p[4][3];
  */
  vtkIdType* pid = new vtkIdType[4];
  double** p = new double*[4];
  for (int i = 0; i < xSize-1; ++i)
    {
      for (int j = 0; j < ySize-1; ++j)
	{
	  // The coordinates of the points in the cell
	  p[0] = new double[3] {x->at(j*xSize+i), y->at(j*xSize+i), z->at(j*xSize+i)};
	  p[1] = new double[3] {x->at(j*xSize+i+1), y->at(j*xSize+i+1), z->at(j*xSize+i+1)};
	  p[2] = new double[3] {x->at((j+1)*xSize+i), y->at((j+1)*xSize+i), z->at((j+1)*xSize+i)};
	  p[3] = new double[3] {x->at((j+1)*xSize+i+1), y->at((j+1)*xSize+i+1), z->at((j+1)*xSize+i+1)};
	  // Insert the points
	  pid[0] = points->InsertNextPoint(p[0]);
	  pid[1] = points->InsertNextPoint(p[1]);
	  pid[2] = points->InsertNextPoint(p[2]);
	  pid[3] = points->InsertNextPoint(p[3]);
	  // Insert the next cell
	  vertices->InsertNextCell(4, pid);
	  // Insert height to colormap
	  colormap->InsertNextValue(z->at(j*xSize +i));
	  colormap->InsertNextValue(z->at(j*xSize +i+1));
	  colormap->InsertNextValue(z->at((j+1)*xSize +i));
	  colormap->InsertNextValue(z->at((j+1)*xSize +i+1));
	  delete[] p[0];
	  p[0] = NULL;
	  delete[] p[1];
	  p[1] = NULL;
	  delete[] p[2];
	  p[2] = NULL;
	  delete[] p[3];
	  p[3] = NULL;
	}
    }
}

void plotSurface(vtkRenderWindowInteractor *iren,
                 std::string filename, int xSize, int ySize)
{
  //int xSize = 0x400;
  //int ySize = 0x400;
  float xMin = -PI;
  float xMax = PI;
  float yMin = -PI;
  float yMax = PI;
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
  /*
  float *x = new float[xSize];
  float *y = new float[ySize];
  float *z = new float[ySize*xSize];
  */
  std::vector<double> *x = new std::vector<double>();
  std::vector<double> *y = new std::vector<double>();
  std::vector<double> *z = new std::vector<double>();
  /*
  createSurfaceArray(x, y, z, xSize, ySize,
                     xMin, xMax, yMin, yMax);
  */
  createSurfaceArrayFromFile(x, y, z, xSize, ySize,
                             xMin, xMax, yMin, yMax,
                             filename);
  /* Find max and min for the colormap */
  std::vector<double> *zBounds = new std::vector<double>();
  zBounds->push_back(z->at(0));
  zBounds->push_back(z->at(0));
  calcZBounds(z, zBounds);

  /* Create the topology of the point (a vertex) */
  createTopology(points, vertices, colormap,
                x, y, z, xSize, ySize);
  /* Clean up unused variables */
  /*
  delete[] x;
  x = NULL;
  delete[] y;
  y = NULL;
  delete[] z;
  z = NULL;
  */
  /*
    ----------------------------------------------------------------------------
    Set the points and vertices created as the geometry and topology of the
    polydata. Create primitive graphics from polygons and addactor it to the
    rendering scene.
    ----------------------------------------------------------------------------
  */
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
  mapper->SetScalarRange(zBounds->at(0), zBounds->at(1)); // colormap boundaries
  //delete[] zBounds;
  //zBounds = NULL;
  /* set the rendering scene */
  actor->SetMapper(mapper);
  actor->GetProperty()->SetOpacity(0.7);
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
  iren->SetRenderWindow(renWin); // set the window for the interactor
}

int main(int argc, char *argv[], char *envp[])
{
  int nxpts, nypts;
  std::string filename;
  if (argc > 3)
    {
      filename = argv[1];
      nxpts = std::stoi(argv[2]);
      nypts = std::stoi(argv[3]);
    }
  std::cout << filename << ", " << nxpts << ", " << nypts << ", " << std::endl;

  for (int i = 0; i < argc; i++)
    {
      std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
    }

  // Interactor to interact with the render window
  VTK_SP(vtkRenderWindowInteractor, iren);
  //std::vector<double> surfData = loadData(filename, nxpts, nypts);
  plotSurface(iren, filename, nxpts, nypts);
  //renWin->Render();
  //ren->GetActiveCamera()->Zoom(1.5);
  iren->Render();
  iren->Start();
  return EXIT_SUCCESS;
}
/*
int functn(float *a, int val)
{
  //std::cout << &a << std::endl;
  for (int i = 0; i < 5; i++)
    {
      //std::cout << a[i] << std::endl;
      a[i] = val++;
    }
}

int functn(std::vector<double> *a, double val)
{
  std::cout << &a << std::endl;
  for (int i = 0; i < 5; i++)
    {
      a->push_back(val);
      val += 1;
    }
  return val;
}

int main(int argc, char *argv[], char *envp[])
{
  double val = 0;
  std::vector<double> *X = new std::vector<double>();
  std::vector<double> *Y = new std::vector<double>();
  std::cout << &x << std::endl;
  val = functn(X, val);
  std::cout << &y << std::endl;
  val = functn(Y, val);
  for (std::vector<double>::iterator i = X->begin(); i != X->end(); i++)
    {
      std::cout << *i << ", ";
    }
  std::cout << "" << std::endl;
  for (std::vector<double>::iterator i = Y->begin(); i != Y->end(); i++)
    {
      std::cout << *i << ", ";
    }
  std::cout << "" << std::endl;
}
*/
