#ifndef SRC_PLOT_HPP_
#define SRC_PLOT_HPP_

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>

#include <string>
#include <vector>
#include <array>

template<typename T>
struct Bounds
{
  T min;
  T max;
};

class Plot
{
public:
  Plot(std::string &fileX, std::string &fileY, std::string &fileZ, int xSize, int ySize)
    : fileX(fileX), fileY(fileY), fileZ(fileZ), xSize(xSize), ySize(ySize)
  {
  }

  vtkSmartPointer<vtkPolyDataMapper> createSurface();
  void plotSurface(vtkSmartPointer<vtkPolyDataMapper> mapper);


private:
  std::string fileX;
  std::string fileY;
  std::string fileZ;
  int xSize;
  int ySize;

  std::vector<float> BytesToDouble(std::vector<char> const &src);
  std::vector<char> ReadAllBytes(std::string const &filename);
  inline std::vector<float> loadData(std::string const &filename)
  {
    return BytesToDouble(ReadAllBytes(filename));
  }
  Bounds<float> calcZBounds(std::vector<float> const &arr);
  Bounds<float> createTopology(vtkSmartPointer<vtkPolyData> point,
			       std::vector<float> const &x,
			       std::vector<float> const &y,
			       std::vector<float> const &z);
  vtkSmartPointer<vtkActor> makeActor(vtkSmartPointer<vtkPolyDataMapper> mapper);
  vtkSmartPointer<vtkRenderer> makeRender(vtkSmartPointer<vtkPolyDataMapper> mapper);
  vtkSmartPointer<vtkRenderWindow> makeRenderWindow(vtkSmartPointer<vtkPolyDataMapper> mapper);
};


#endif /* SRC_PLOT_HPP_ */
