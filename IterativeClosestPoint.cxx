/*=========================================================================
  Source:    Examples/Registration/IterativeClosestPoint1.cxx
=========================================================================*/
#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

// ITK Includes
#include <itkCenteredRigid2DTransform.h>
#include <itkEuclideanDistancePointMetric.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include <itkPointSet.h>
#include <itkPointSetToPointSetRegistrationMethod.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkImageToVTKImageFilter.h"

// VTK Includes
#include <vtkSmartPointer.h>
#include <vtkActor2D.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageViewer.h>
#include <vtkPointPicker.h>
#include <vtkCommand.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPoints.h>
#include <vtkPointSource.h>
#include <vtkPointSetSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkProperty2D.h>
#include <vtkGlyph2D.h>
#include <vtkOutlineSource.h>
#include <vtkTriangleStrip.h>
#include <vtkCellArray.h>

// FLTK Includes
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>

#include <iostream>
#include <fstream>
#include <limits.h>

vtkSmartPointer<vtkPointSource> pointSource;

class vtkButtonPressCallback : public vtkCommand
{
private:
	vtkPoints *pointsList;
	vtkPolyData *polyData;
	vtkRenderWindow *renderWindow;
	bool inPolygon;

public: 
	vtkButtonPressCallback() {
		pointsList = vtkPoints::New();
		polyData = vtkPolyData::New();
		polyData->Allocate();
		inPolygon = false;
	}
	static vtkButtonPressCallback *New() { return new vtkButtonPressCallback; }

	virtual vtkPolyData *GetPolyData() {
		return polyData;
	}

	virtual void SetRenderWindow(vtkRenderWindow *window) {
		this->renderWindow = window;
	}
	
	virtual void HandleLeftClick(vtkRenderWindowInteractor *interactor) {
		if(!inPolygon) {
			polyData->Allocate();
			pointsList->Delete();
			pointsList = vtkPoints::New();
		}
		pointsList->InsertNextPoint(interactor->GetEventPosition()[0], 
			interactor->GetEventPosition()[1], 
			0);
		polyData->SetPoints(pointsList);

		std::cout << "( " << interactor->GetEventPosition()[0] <<
        	         ", " << interactor->GetEventPosition()[1] << ");";

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

		inPolygon = true;
	}

	virtual void HandleRightClick(vtkRenderWindowInteractor *interactor) {
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
			polyData->InsertNextCell(VTK_TRIANGLE_STRIP, numPoints, polyConnection);
			vtkTriangleStrip *cell = 
				reinterpret_cast<vtkTriangleStrip*>(polyData->GetCell(polyData->GetNumberOfCells() - 1));

			polyData->Print(std::cout);
			std::cout << " Triangle Strips: " << std::endl;
			polyData->GetStrips()->Print(std::cout);

			std::cout << " length(cells) = "
				  << polyData->GetNumberOfCells() << std::endl;

			delete []polyConnection;
		}

		inPolygon = false;
	}

	virtual void Execute(vtkObject *caller, unsigned long eventType, void *)
	{
		vtkRenderWindowInteractor *interactor = 
			reinterpret_cast<vtkRenderWindowInteractor*>(caller);

		if(eventType == LeftButtonPressEvent)
			HandleLeftClick(interactor);
		else if(eventType == RightButtonPressEvent)
			HandleRightClick(interactor);
	}
};

class CommandIterationUpdate : public itk::Command
{
public:
	typedef CommandIterationUpdate Self;
	typedef itk::Command Superclass;
	typedef itk::SmartPointer <Self> Pointer;
	itkNewMacro(Self);

protected:
	CommandIterationUpdate ()
	{
	};

public:
	typedef itk::LevenbergMarquardtOptimizer OptimizerType;
	typedef const OptimizerType *OptimizerPointer;
	
	void Execute(itk::Object * caller, const itk::EventObject & event)
	{
		Execute((const itk::Object *) caller, event);
	}

	void Execute(const itk::Object * object, const itk::EventObject & event)
  	{
    	OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>(object);

		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
      	}

		std::cout << "Difference = " << optimizer->GetCachedValue() << std::endl;
    	std::cout << "Position   = " << optimizer->GetCachedCurrentPosition();
    	std::cout << std::endl;
	}
};


int main (int argc, char *argv[])
{
	if (argc < 6)
	{
    	std::cerr << "Arguments Missing. " << std::endl;
    	std::cerr << "Usage:  IterativeClosestPoint1   fixedPointsFile  movingPointsFile  "
				  << "inputFixedImage  inputMovingImage  outputImage"
				  << std::endl;
      	return 1;
    }

	const unsigned int Dimension = 2;

	typedef itk::PointSet < float, Dimension > PointSetType;

	PointSetType::Pointer fixedPointSet = PointSetType::New();
	PointSetType::Pointer movingPointSet = PointSetType::New();

	typedef PointSetType::PointType PointType;
	typedef PointSetType::PointsContainer PointsContainer;
	typedef PointsContainer::Iterator PointsIterator;

	typedef unsigned char PixelType;
	typedef itk::Image< PixelType, Dimension > InputImageType;
	typedef itk::ImageFileReader< InputImageType > ReaderType;
	typedef itk::Image< PixelType, Dimension > OutputImageType;
	typedef itk::ImageFileWriter< OutputImageType> WriterType;
	typedef itk::ExtractImageFilter<InputImageType, InputImageType> ExtractFilterType;
	typedef itk::ResampleImageFilter<InputImageType, InputImageType> ResampleFilterType;
	typedef itk::SubtractImageFilter<InputImageType, InputImageType, InputImageType> SubtractFilterType;
	typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, double> InterpolatorType;
	typedef itk::ImageToVTKImageFilter<InputImageType> VTKExporterType;

	PointsContainer::Pointer fixedPointContainer = PointsContainer::New();
	PointsContainer::Pointer movingPointContainer = PointsContainer::New();

	PointType fixedPoint;
	PointType movingPoint;

	// Read the file containing coordinates of fixed points.
	std::ifstream fixedFile;
	fixedFile.open(argv[1]);
  	if(fixedFile.fail())
    {
		std::cerr << "Error opening points file with name : " << std::endl;
		std::cerr << argv[1] << std::endl;
		return 2;
    }

	unsigned int pointId = 0;
	fixedFile >> fixedPoint;
	while (!fixedFile.eof())
	{
    	fixedPointContainer->InsertElement(pointId, fixedPoint);
    	fixedFile >> fixedPoint;
    	pointId++;
    }
	fixedPointSet->SetPoints(fixedPointContainer);
	std::cout << "Number of fixed Points = " <<
    	fixedPointSet->GetNumberOfPoints() << std::endl;

	// Read the file containing coordinates of moving points.
	std::ifstream movingFile;
	movingFile.open (argv[2]);
	if (movingFile.fail())
    {
    	std::cerr << "Error opening points file with name : " << std::endl;
    	std::cerr << argv[2] << std::endl;
    	return 2;
    }

	pointId = 0;
	movingFile >> movingPoint;
	while (!movingFile.eof())
    {
	    movingPointContainer->InsertElement(pointId, movingPoint);
    	movingFile >> movingPoint;
	    pointId++;
    }
	movingPointSet->SetPoints(movingPointContainer);
	std::cout << "Number of moving Points = "
    	<< movingPointSet->GetNumberOfPoints() << std::endl;

	//------------------------------------------------------------
	// Read input images
	//------------------------------------------------------------
	ReaderType::Pointer fixedReader = ReaderType::New();
	ReaderType::Pointer movingReader = ReaderType::New();
	std::cout << "Reading fixed image from " << argv[3] << std::endl;
	std::cout << "Reading moving image from " << argv[4] << std::endl;
	fixedReader->SetFileName(argv[3]);
	movingReader->SetFileName(argv[4]);
	try {
		fixedReader->Update();
		movingReader->Update();
	}
	catch(itk::ExceptionObject &ex) {
    	std::cout << ex << std::endl;
	    return EXIT_FAILURE;
    }

	//-----------------------------------------------------------
	// Create Bounding Box
	//-----------------------------------------------------------
	PointsIterator pointIterator = fixedPointContainer->Begin();
	PointsIterator end = fixedPointContainer->End();
	
	float fixedXMin = INT_MAX;
	float fixedYMin = INT_MAX;
	float fixedXMax = 0;
	float fixedYMax = 0;

	while( pointIterator != end )
	{
		PointType p = pointIterator.Value();
		if(p[0] < fixedXMin) fixedXMin = p[0];
		if(p[0] > fixedXMax) fixedXMax = p[0];
		if(p[1] < fixedYMin) fixedYMin = p[1];
		if(p[1] > fixedYMax) fixedYMax = p[1];
		++pointIterator;
	}
	std::cout << "Fixed Bounding Box UL:(" << fixedXMin << ", " << fixedYMin <<
			  ") BR:(" << fixedXMax << ", " << fixedYMax << ")" << std::endl;
	ExtractFilterType::Pointer extractor = ExtractFilterType::New();
	ExtractFilterType::InputImageRegionType region;
	region.SetIndex(0, fixedXMin);
	region.SetSize(0, fixedXMax - fixedXMin);
	region.SetIndex(1, fixedYMin);
	region.SetSize(1, fixedYMax - fixedYMin);
	std::cout << region << std::endl;
	extractor->SetExtractionRegion(region);

	//-----------------------------------------------------------
	// Set up Image Writer
	//-----------------------------------------------------------
	WriterType::Pointer resultWriter = WriterType::New();
	resultWriter->SetFileName(argv[5]);
	std::cout << "Writing output image to " << argv[5] << std::endl;
	//resultWriter->SetInput(extractor->GetOutput());
	//resultWriter->Update();

	//-----------------------------------------------------------
	// Set up  the Metric
	//-----------------------------------------------------------
	typedef itk::EuclideanDistancePointMetric<PointSetType, PointSetType>
	    MetricType;
	typedef MetricType::TransformType TransformBaseType;
  	typedef TransformBaseType::ParametersType ParametersType;
	typedef TransformBaseType::JacobianType JacobianType;

	MetricType::Pointer metric = MetricType::New();

	//-----------------------------------------------------------
	// Set up a Transform
	//-----------------------------------------------------------
	typedef itk::CenteredRigid2DTransform < double > TransformType;
	TransformType::Pointer transform = TransformType::New();

	// Optimizer Type
	typedef itk::LevenbergMarquardtOptimizer OptimizerType;

	OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetUseCostFunctionGradient(false);

	// Registration Method
	typedef itk::PointSetToPointSetRegistrationMethod 
		<PointSetType, PointSetType> RegistrationType;

	RegistrationType::Pointer registration = RegistrationType::New();

	// Scale the translation components of the Transform in the Optimizer
	OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
	scales[0] = itk::Math::pi/1000.0;
	scales[1] = 1.0/10.0;
	scales[2] = 1.0/10.0;
	scales[3] = 1.0/100.0;
	scales[4] = 1.0/100.0;
	
	ParametersType defaults(transform->GetNumberOfParameters());

	double guess_angle = 3*itk::Math::pi_over_2;
	PointType v1_p1 = fixedPointContainer->GetElement(0);
	PointType v1_p2 = fixedPointContainer->GetElement(1);
	double v1_x = v1_p2[0] - v1_p1[0];
	double v1_y = v1_p2[1] - v1_p1[1];
	PointType v2_p1 = movingPointContainer->GetElement(0);
	PointType v2_p2 = movingPointContainer->GetElement(1);
	double v2_x = v2_p2[0] - v2_p1[0];
	double v2_y = v2_p2[1] - v2_p1[1];
	guess_angle = atan2(v1_y,v1_x) - atan2(v2_y,v2_x);
	std::cout << "guess_angle = " << guess_angle << std::endl;

	// rotation
	defaults.SetElement(0, guess_angle);
	// center x,y
	defaults.SetElement(1, 0.0);
	defaults.SetElement(2, 0.0);
	// translation x,y 
	defaults.SetElement(3, 0.0);
	defaults.SetElement(4, 0.0);
	transform->SetParameters(defaults);

	unsigned long numberOfIterations = 500;
	double gradientTolerance = 1e-5;	// convergence criterion
	double valueTolerance = 1e-5;	// convergence criterion
	double epsilonFunction = 1e-6;	// convergence criterion

	optimizer->SetScales(scales);
	optimizer->SetNumberOfIterations(numberOfIterations);
	optimizer->SetValueTolerance(valueTolerance);
	optimizer->SetGradientTolerance(gradientTolerance);
	optimizer->SetEpsilonFunction(epsilonFunction);

	registration->SetInitialTransformParameters(transform->GetParameters());

	//------------------------------------------------------
	// Connect all the components required for Registration
	//------------------------------------------------------
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetTransform(transform);
	registration->SetFixedPointSet(fixedPointSet);
	registration->SetMovingPointSet(movingPointSet);

	// Connect an observer
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	try
	{
		registration->StartRegistration();
	}
  	catch (itk::ExceptionObject & e)
	{
    	std::cout << e << std::endl;
	    return EXIT_FAILURE;
	}

	std::cout << "Solution = " << transform->GetParameters() << std::endl;

	//------------------------------------------------------
	// Transform and extract solution from moving image
	//------------------------------------------------------
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	InputImageType::ConstPointer fixedImage = fixedReader->GetOutput();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	InputImageType::SizeType outputImageSize;
	outputImageSize[0] = 500;
	outputImageSize[1] = 500;
	
	TransformType::Pointer inverseTransform = TransformType::New();
	transform->GetInverse(inverseTransform);

	resampler->SetSize(outputImageSize);
	resampler->SetTransform(inverseTransform);
	resampler->SetDefaultPixelValue(0);
	resampler->SetInterpolator(interpolator);
	resampler->SetOutputOrigin(fixedImage->GetOrigin());
	resampler->SetOutputSpacing(fixedImage->GetSpacing());
	resampler->SetOutputDirection(fixedImage->GetDirection());
	
	SubtractFilterType::Pointer subtractor = SubtractFilterType::New();	
	subtractor->SetInput1(fixedImage);
	
	resampler->SetInput(movingReader->GetOutput());
	subtractor->SetInput2(resampler->GetOutput());
	extractor->SetInput(subtractor->GetOutput());
	resultWriter->SetInput(extractor->GetOutput());

	//------------------------------------------------------
	// Export image to VTK and display
	//------------------------------------------------------
	VTKExporterType::Pointer exporter = VTKExporterType::New();
	exporter->SetInput(extractor->GetOutput());
	
	vtkSmartPointer<vtkImageViewer> viewer = vtkImageViewer::New();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkRenderWindowInteractor::New();
	vtkSmartPointer<vtkPointPicker> pointPicker = vtkPointPicker::New();
	vtkSmartPointer<vtkButtonPressCallback> callback = vtkButtonPressCallback::New();
	vtkSmartPointer<vtkGlyph2D> glyph = vtkGlyph2D::New();
	vtkSmartPointer<vtkPolyDataMapper2D> pointMapper = vtkPolyDataMapper2D::New();
	vtkSmartPointer<vtkPolyDataMapper2D> lineMapper = vtkPolyDataMapper2D::New();
	vtkSmartPointer<vtkActor2D> coordinatesActor = vtkActor2D::New();
	vtkSmartPointer<vtkActor2D> linesActor = vtkActor2D::New();
	vtkSmartPointer<vtkProperty2D> coordinatesProp = coordinatesActor->GetProperty();
	vtkSmartPointer<vtkProperty2D> linesProp = linesActor->GetProperty();
	pointSource = vtkPointSource::New();

	coordinatesProp->SetColor(1.0, 0.0, 0.0);
	coordinatesProp->SetPointSize(3.0);
	linesProp->SetColor(0.8, 0.0, 0.0);
	linesProp->SetOpacity(0.3);

	viewer->SetupInteractor(renderWindowInteractor);
	renderWindowInteractor->SetPicker(pointPicker);
	callback->SetRenderWindow(viewer->GetRenderWindow());
  
	glyph->SetSource(pointSource->GetOutput());
	glyph->SetInput(callback->GetPolyData());
	pointMapper->SetInputConnection(glyph->GetOutputPort());
	coordinatesActor->SetMapper(pointMapper);

	lineMapper->SetInput(callback->GetPolyData());
	linesActor->SetMapper(lineMapper);

	vtkSmartPointer<vtkRenderer> renderer = viewer->GetRenderer();
	renderer->AddActor(coordinatesActor);
	renderer->AddActor(linesActor);
	renderer->SetBackground(0.8, 0.8, 0.8);
	renderWindowInteractor->AddObserver(
		vtkCommand::LeftButtonPressEvent, 
		callback
	);
	renderWindowInteractor->AddObserver(
		vtkCommand::RightButtonPressEvent, 
		callback
	);

	viewer->SetInput(exporter->GetOutput());
	viewer->SetSize(500,500);
	viewer->SetColorWindow(255);
	viewer->SetColorLevel(128);
	
	viewer->Render();
	renderWindowInteractor->Start();
		
	try
	{
		resultWriter->Update();
	}
  	catch (itk::ExceptionObject & e)
	{
    	std::cout << e << std::endl;
	    return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
