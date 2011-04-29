#include <vtkPoints.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

class vtkPointSelectorStyle : public vtkInteractorStyleTrackballCamera
{
private:
	vtkSmartPointer<vtkPoints> pointsList;
	vtkSmartPointer<vtkPolyData> polyData;
	vtkSmartPointer<vtkPointSource> pointSource;
	bool inPolygon;

public:
	vtkPointSelectorStyle();
	static vtkPointSelectorStyle *New();
	virtual vtkPolyData *GetPolyData();
	virtual void SetPointSource(vtkPointSource *pointSource);
	virtual void OnLeftButtonDown();
	virtual void OnRightButtonDown();
	virtual vtkPoints *GetPoints();
};
