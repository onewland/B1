#include "vtkItkSourceViewer.h"

template <class TInputImage>
vtkITKSourceViewer<TInputImage>::vtkITKSourceViewer(
	int x, int y, int w, int h) : vtkFlRenderWindowInteractor(x,y,w,h,"")
{
	viewer = vtkImageViewer::New();
	renderer = viewer->GetRenderer();
	exporter = VTKExporterType::New();
	callback = vtkPointSelectorStyle::New();
	glyph = vtkGlyph2D::New();
	pointMapper = vtkPolyDataMapper2D::New();
	lineMapper = vtkPolyDataMapper2D::New();
	coordinatesActor = vtkActor2D::New();
	linesActor = vtkActor2D::New();
	pointSource = vtkPointSource::New();
	verticalFlipper = FlipFilterType::New();
	windowingFilter = WindowFilterType::New();
	brightnessFilter = AddConstantToImageFilterType::New();
	
	// Flip across Y, but not X
	bool axes[] = {false, true};
	verticalFlipper->SetFlipAxes(axes);

	coordinatesProp = coordinatesActor->GetProperty();
	linesProp = linesActor->GetProperty();

	callback->SetPointSource(pointSource);
	viewer->SetupInteractor(this);
	this->SetInteractorStyle(callback);

	SetRenderWindow(viewer->GetRenderWindow());
	
	windowingFilter->SetInput(verticalFlipper->GetOutput());
	brightnessFilter->SetInput(windowingFilter->GetOutput());
	exporter->SetInput(brightnessFilter->GetOutput());

	viewer->SetInput(exporter->GetOutput());
	viewer->SetColorWindow(255);
	viewer->SetColorLevel(128);

	coordinatesProp->SetColor(1.0, 0.0, 0.0);
	coordinatesProp->SetPointSize(3.0);
	linesProp->SetColor(0.8, 0.0, 0.0);
	linesProp->SetOpacity(0.3);

	glyph->SetSource(pointSource->GetOutput());
	glyph->SetInput(callback->GetPolyData());
	pointMapper->SetInputConnection(glyph->GetOutputPort());
	coordinatesActor->SetMapper(pointMapper);
	
	lineMapper->SetInput(callback->GetPolyData());
	linesActor->SetMapper(lineMapper);

	renderer->AddActor(coordinatesActor);
	renderer->AddActor(linesActor);
	renderer->SetBackground(0.8, 0.8, 0.8);
}

template <class TInputImage>
void vtkITKSourceViewer<TInputImage>::SetImageInput(const TInputImage *input)
{
	verticalFlipper->SetInput(input);
}

template <class TInputImage>
vtkImageViewer* vtkITKSourceViewer<TInputImage>::GetViewer()
{
	return viewer;
}

template <class TInputImage>
vtkPoints* vtkITKSourceViewer<TInputImage>::GetSelectedPoints()
{
	return callback->GetPoints();
}

template <class TInputImage>
InputImageFiltersDialog<TInputImage> *
vtkITKSourceViewer<TInputImage>::GetFilterDialog()
{
	return filterDialog;
}

template <class TInputImage>
void
vtkITKSourceViewer<TInputImage>::MakeFilterDialog(const char *title = "Pre-processing options")
{
	filterDialog = new FilterDialogType(title);
	filterDialog->SetInputViewer(this);
}

template <class TInputImage>
void 
vtkITKSourceViewer<TInputImage>::SetIntensityWindowMin
(double intensity)
{
	windowingFilter->SetWindowMinimum((unsigned short) intensity);
	UpdateView();
}

template <class TInputImage>
void 
vtkITKSourceViewer<TInputImage>::SetIntensityWindowMax
(double intensity)
{
	windowingFilter->SetWindowMaximum((unsigned short) intensity);
	UpdateView();
}

template <class TInputImage>
void 
vtkITKSourceViewer<TInputImage>::SetIntensityOffset
(double intensity)
{
	brightnessFilter->GetFunctor().SetOffset(intensity);
	brightnessFilter->Modified();
	UpdateView();
}

template <class TInputImage>
void
vtkITKSourceViewer<TInputImage>::UpdateView()
{
	viewer->Render();
}
