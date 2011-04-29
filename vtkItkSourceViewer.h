#ifndef __vtkItkSourceViewer_h
#define __vtkItkSourceViewer_h

#include "itkImageToVTKImageFilter.h"
#include "vtkFlRenderWindowInteractor.h"
#include "vtkPointSelectorStyle.h"
#include "InputImageFiltersDialog.h"

#include <itkImageSource.h>
#include <itkFlipImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>

#include <vtkSmartPointer.h>
#include <vtkImageViewer.h>
#include <vtkRenderer.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>
#include <vtkGlyph2D.h>
#include <vtkActor2D.h>

template <class TInput>
class AddConstantFunctor
{
private:
	typedef TInput IntensityVal;
	IntensityVal offset;

public:
	AddConstantFunctor() {
		offset = 0;
	}
	
	void SetOffset(IntensityVal newOffset) {
		offset = newOffset;
	}

	inline TInput operator()(const TInput &in) const {
		return (in + offset > 255 ? 255 : in + offset);
	}
};

template <class TInputImage> 
class vtkITKSourceViewer : public vtkFlRenderWindowInteractor {
public:
	typedef TInputImage InputImageType;
	typedef typename InputImageType::PixelType PixelType;
	typedef itk::ImageToVTKImageFilter<InputImageType> VTKExporterType;
	typedef itk::FlipImageFilter<InputImageType> FlipFilterType;
	typedef itk::IntensityWindowingImageFilter<InputImageType> WindowFilterType;
	typedef InputImageFiltersDialog<InputImageType> FilterDialogType;
	typedef AddConstantFunctor<PixelType> AddConstantFunctorType;
	typedef itk::UnaryFunctorImageFilter<InputImageType, InputImageType, AddConstantFunctorType>
		AddConstantToImageFilterType;

private:
	vtkSmartPointer<vtkImageViewer> viewer;
	vtkSmartPointer<vtkRenderer> renderer;
	typename VTKExporterType::Pointer exporter;
	typename FlipFilterType::Pointer verticalFlipper;
	typename WindowFilterType::Pointer windowingFilter;
	typename AddConstantToImageFilterType::Pointer brightnessFilter;
	vtkSmartPointer<vtkPointSelectorStyle> callback;
	vtkSmartPointer<vtkGlyph2D> glyph;
	vtkSmartPointer<vtkPolyDataMapper2D> pointMapper;
	vtkSmartPointer<vtkPolyDataMapper2D> lineMapper;
	vtkSmartPointer<vtkActor2D> coordinatesActor;
	vtkSmartPointer<vtkActor2D> linesActor;
	vtkSmartPointer<vtkProperty2D> coordinatesProp;
	vtkSmartPointer<vtkProperty2D> linesProp;
	vtkSmartPointer<vtkPointSource> pointSource;
	FilterDialogType *filterDialog;

public:
	vtkITKSourceViewer(int x, int y, int w, int h);

	virtual void SetImageInput(const InputImageType *input);
	virtual vtkImageViewer *GetViewer();
	virtual vtkPoints *GetSelectedPoints();
	virtual FilterDialogType *GetFilterDialog();
	virtual void MakeFilterDialog(const char *title);
	virtual void SetIntensityWindowMin(double intensity);
	virtual void SetIntensityWindowMax(double intensity);
	virtual void SetIntensityOffset(double offset);
	virtual void UpdateView();
};

#endif
