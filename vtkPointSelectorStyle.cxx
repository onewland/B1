#include "vtkPointSelectorStyle.h"

vtkPointSelectorStyle::vtkPointSelectorStyle() {
	pointsList = vtkPoints::New();
	polyData = vtkPolyData::New();
	polyData->Allocate();
	inPolygon = false;
}

vtkPointSelectorStyle *vtkPointSelectorStyle::New() { 
	return new vtkPointSelectorStyle; 
}

vtkPolyData *vtkPointSelectorStyle::GetPolyData() {
	return polyData;
}

void vtkPointSelectorStyle::OnLeftButtonDown() {
	vtkSmartPointer<vtkRenderWindowInteractor> interactor = GetInteractor();
	if(!inPolygon) {
		polyData->Allocate();
		pointsList->Delete();
		pointsList = vtkPoints::New();
	}
	std::cout << "( " << interactor->GetEventPosition()[0] <<
       	         ", " << interactor->GetEventPosition()[1] << ");";

	pointsList->InsertNextPoint(interactor->GetEventPosition()[0], 
		interactor->GetEventPosition()[1], 
		0);
	polyData->SetPoints(pointsList);

	int lastPointIndex = pointsList->GetNumberOfPoints() - 1;
	if(inPolygon) {
		vtkIdType connection[2];
		connection[0] = lastPointIndex;
		connection[1] = lastPointIndex - 1;
		polyData->InsertNextCell(VTK_LINE, 2, connection);
	}

	vtkIdType pointIndex[1];
	pointIndex[0] = lastPointIndex;
	pointSource->SetNumberOfPoints(lastPointIndex + 1);
	
	std::cout << " length(pointsList) = " 
			  << pointsList->GetNumberOfPoints() << std::endl;

	std::cout << " length(cells) = "
			  << polyData->GetNumberOfCells() << std::endl;

	pointSource->Update();

	inPolygon = true;
	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void vtkPointSelectorStyle::OnRightButtonDown() {
	vtkSmartPointer<vtkRenderWindowInteractor> interactor = GetInteractor();
	if(inPolygon) {
		int numPoints = pointsList->GetNumberOfPoints();
		vtkIdType connection[2];
		vtkIdType *polyConnection = new vtkIdType[numPoints];
		for(int i = 0; i < numPoints; i++) {
			polyConnection[i] = i;
		}
		connection[0] = pointsList->GetNumberOfPoints() - 1;
		connection[1] = 0;
		polyData->InsertNextCell(VTK_LINE, 2, connection);

		// Note that I was trying to use VTK_TRIANGLE_STRIP to fill the interior of
		// the region created by point selection. VTK_POLY doesn't work for 
		// irregular polygons which are likely to appear. To fill in the
		// selection would probably require a custom region filling PolyDataMapper
		// but VTK_TRIANGLE_STRIP doesn't work (like this) -- I just didn't delete
		// the code because ostensibly any PolyDataMapper will have to have all
		// the points in a single cell to do the region fill.
		polyData->InsertNextCell(VTK_TRIANGLE_STRIP, numPoints, polyConnection);

		std::cout << " length(cells) = "
			  << polyData->GetNumberOfCells() << std::endl;

		delete []polyConnection;
	}

	inPolygon = false;
	vtkInteractorStyleTrackballCamera::OnRightButtonDown();
}

void vtkPointSelectorStyle::SetPointSource(vtkPointSource *source) {
	pointSource = source;
}

vtkPoints *vtkPointSelectorStyle::GetPoints() {
	return pointsList;
}
