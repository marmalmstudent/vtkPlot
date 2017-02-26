/*
  This code is optimized for 80 columns
|------------------------------------------------------------------------------|
 */
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

static void BytesToDouble(std::vector<char> *src, std::vector<double> *dst,
                          bool littleEndian)
{
  /*
    ----------------------------------------------------------------------------
    This function loads an array of bytes and assembles it into an array of
    doubles.
    ----------------------------------------------------------------------------
   */
  int idx = 0;
  int const doubleSize = sizeof(double);
  int srcSize = src->size();
  double val;
  union {
    double myDouble;
    char myChars[doubleSize];
  } uni;
  if (littleEndian) {
    while(idx < srcSize) {
      for(int k = 0; k < doubleSize; k++)
	uni.myChars[k] = src->at(idx+k);
      dst->push_back(uni.myDouble);
      idx += doubleSize;
    }
  }
  else {
    while(idx < srcSize) {
      for(int k = 0; k < doubleSize; k++)
	uni.myChars[k] = src->at(idx + doubleSize-1 - k);
      dst->push_back(uni.myDouble);
      idx += doubleSize;
    }
  }
}

static void ReadAllBytes(char const* filename, std::vector<char>  *result)
{
  /*
    ----------------------------------------------------------------------------
    This function reads bytes from a file and stores the data in the
    supplied byte array.
    ----------------------------------------------------------------------------
   */
  std::ifstream ifs(filename, std::ifstream::binary|std::ifstream::ate);
  if (ifs.is_open()) {
    std::ifstream::pos_type pos = ifs.tellg();
    
    result->resize(pos);
    
    ifs.seekg(0, std::ios::beg);
    ifs.read(&(*result)[0], pos);  // the adress of result[0]
    ifs.close();
  }
  else {
    std::cout << "unable to find file: " << filename << std::endl;
  }
  return;
}

void loadData(std::string filename, std::vector<double> *arr)
{
  /*
    ----------------------------------------------------------------------------
    This function reads bytes from a file and assembles it into an array
    of doubles.
    ----------------------------------------------------------------------------
   */
  std::vector<char> b;
  ReadAllBytes(filename.c_str(), &b);
  BytesToDouble(&b, arr, true);
}

void calcZBounds(std::vector<double> *arr, std::vector<double> *zBounds)
{
  /*
    ----------------------------------------------------------------------------
    Find then maximum and minimum values of a given array to calibrate colormap.
    ----------------------------------------------------------------------------
   */
  for (std::vector<double>::iterator i = arr->begin(); i != arr->end(); i++) {
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
  vtkIdType *pid = new vtkIdType[4];
  for (int i = 0; i < xSize-1; ++i) {
      for (int j = 0; j < ySize-1; ++j) {
	// The coordinates of the points in the cell
	// Insert the points
	pid[0] = points->InsertNextPoint(x->at(j*xSize+i),
					 y->at(j*xSize+i),
					 z->at(j*xSize+i));
	pid[1] = points->InsertNextPoint(x->at(j*xSize+i+1),
					 y->at(j*xSize+i+1),
					 z->at(j*xSize+i+1));
	pid[2] = points->InsertNextPoint(x->at((j+1)*xSize+i),
					 y->at((j+1)*xSize+i),
					 z->at((j+1)*xSize+i));
	pid[3] = points->InsertNextPoint(x->at((j+1)*xSize+i+1),
					 y->at((j+1)*xSize+i+1),
					 z->at((j+1)*xSize+i+1));
	// Insert the next cell
	vertices->InsertNextCell(4, pid);
	// Insert height to colormap
	colormap->InsertNextValue(z->at(j*xSize +i));
	colormap->InsertNextValue(z->at(j*xSize +i+1));
	colormap->InsertNextValue(z->at((j+1)*xSize +i));
	colormap->InsertNextValue(z->at((j+1)*xSize +i+1));
      }
  }
  delete[] pid;
  pid = NULL;
}

void plotSurface(vtkRenderWindowInteractor *iren,
                 std::string fileX, std::string fileY,
                 std::string fileZ, int xSize, int ySize)
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
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  loadData(fileX, &x);
  loadData(fileY, &y);
  loadData(fileZ, &z);
  /* Find max and min for the colormap */
  std::vector<double> zBounds;
  zBounds.push_back(z.at(0));
  zBounds.push_back(z.at(0));
  calcZBounds(&z, &zBounds);

  /* Create the topology of the point (a vertex) */
  createTopology(points, vertices, colormap,
                &x, &y, &z, xSize, ySize);
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
  mapper->SetScalarRange(zBounds.at(0), zBounds.at(1)); // colormap boundaries
  //delete[] zBounds;
  //zBounds = NULL;
  /* set the rendering scene */
  actor->SetMapper(mapper);
  actor->GetProperty()->SetOpacity(0.3);
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
  std::string fileX, fileY, fileZ;
  if (argc > 5)
    {
      fileX = argv[1];
      fileY = argv[2];
      fileZ = argv[3];
      nxpts = std::stoi(argv[4]);
      nypts = std::stoi(argv[5]);
    }
  std::cout << "Reading from file: \n\t" << fileX << "\n\t" << fileY << "\n\t" << fileZ << std::endl;
  std::cout << "Grid size: " << nxpts << " x " << nypts << std::endl;
  // Interactor to interact with the render window
  VTK_SP(vtkRenderWindowInteractor, iren);
  //std::vector<double> surfData = loadData(filename, nxpts, nypts);
  plotSurface(iren, fileX, fileY, fileZ, nxpts, nypts);
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
