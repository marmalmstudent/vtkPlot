#ifndef SRC_PLOT_HPP_
#define SRC_PLOT_HPP_

#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderer.h>
#include <vtk/vtkActor.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPolyDataMapper.h>

#include <string>
#include <vector>
#include <array>

class Plot
{
public:
	Plot(std::string& fileX, std::string& fileY, std::string& fileZ, int xSize, int ySize) :
		fileX(fileX), fileY(fileY), fileZ(fileZ), xSize(xSize), ySize(ySize) { }

	vtkSmartPointer<vtkPolyDataMapper> createSurface();
	void plotSurface(vtkSmartPointer<vtkPolyDataMapper> mapper);


private:
	std::string fileX;
	std::string fileY;
	std::string fileZ;
	int xSize;
	int ySize;

	std::vector<float> BytesToDouble(const std::vector<char> &src);
	std::vector<char> ReadAllBytes(const std::string &filename);
	inline std::vector<float> loadData(const std::string &filename) {
		return BytesToDouble(ReadAllBytes(filename));
	}
	std::array<float, 2> calcZBounds(const std::vector<float> &arr);
	std::array<float, 2> createTopology(vtkSmartPointer<vtkPolyData> point, const std::vector<float> &x, const std::vector<float> &y, const std::vector<float> &z);
	vtkSmartPointer<vtkActor> makeActor(vtkSmartPointer<vtkPolyDataMapper> mapper);
	vtkSmartPointer<vtkRenderer> makeRender(vtkSmartPointer<vtkPolyDataMapper> mapper);
	vtkSmartPointer<vtkRenderWindow> makeRenderWindow(vtkSmartPointer<vtkPolyDataMapper> mapper);
};


#endif /* SRC_PLOT_HPP_ */
