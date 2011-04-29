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
#include <itkSubtractImageFilter.h>
#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkImageToVTKImageFilter.h"
#include "itkRegistrationPipelineManager.h"

// VTK Includes
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageViewer.h>
#include <vtkPoints.h>
#include "vtkItkSourceViewer.cxx"

// FLTK Includes
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>

#include "vtkFlRenderWindowInteractor.h"
#include "Application.h"

#include <iostream>
#include <fstream>
#include <limits.h>

typedef itk::CenteredRigid2DTransform< double > TransformType;
//------------------------------------------------------
// Transform and extract solution from moving image
//------------------------------------------------------
template<class TImageType, class TImageFilter, class TCombineFilter>
void 
CombineImages(const TImageType *fixedImage, 
			   const TImageType *movingImage,
			   const TransformType *transform,
			   TImageFilter *outputImageFilter)
{
	typedef typename itk::ResampleImageFilter<TImageType, TImageType> ResampleFilterType;
	typedef typename ResampleFilterType::Pointer ResampleFilterTypePointer;
	typedef typename itk::NearestNeighborInterpolateImageFunction<TImageType, double> InterpolatorType;
	typedef typename InterpolatorType::Pointer InterpolatorTypePointer;
	typedef typename TImageType::SizeType ImageSizeType;
	typedef TCombineFilter CombineFilterType;
	typedef typename CombineFilterType::Pointer CombineFilterTypePointer;

	ResampleFilterTypePointer resampler = ResampleFilterType::New();
	InterpolatorTypePointer interpolator = InterpolatorType::New();
	CombineFilterTypePointer combiner = CombineFilterType::New();	

	ImageSizeType outputImageSize;
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
	

	combiner->SetInput1(fixedImage);
	resampler->SetInput(movingImage);
	combiner->SetInput2(resampler->GetOutput());
	outputImageFilter->SetInput(combiner->GetOutput());
    outputImageFilter->Update();
}

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
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType> WriterType;
	typedef itk::RegistrationPipelineManager<ImageType, Dimension> RegistrationManager;
	typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractFilterType;
	typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType>
		AbsDiffFilterType;

	RegistrationManager::Pointer registrationManager = RegistrationManager::New();

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
	// Set up Image Writer
	//-----------------------------------------------------------
	WriterType::Pointer resultWriter = WriterType::New();
	resultWriter->SetFileName(argv[5]);
	std::cout << "Writing output image to " << argv[5] << std::endl;

	try
	{
		registrationManager->SetRegistrationPoints(fixedPointSet, movingPointSet);
		registrationManager->PerformRegistration();
	}
  	catch (itk::ExceptionObject & e)
	{
    	std::cout << e << std::endl;
	    return EXIT_FAILURE;
	}

	std::cout << "Solution = " << registrationManager->GetTransform()->GetParameters() << std::endl;

	//------------------------------------------------------
	// Transform and extract solution from moving image
	//------------------------------------------------------
	ImageType::ConstPointer fixedImage = fixedReader->GetOutput();
	ImageType::ConstPointer movingImage = movingReader->GetOutput();

	ImageType::SizeType outputImageSize;
	outputImageSize[0] = 500;
	outputImageSize[1] = 500;
	
	TransformType::Pointer transform = registrationManager->GetTransform();

	//------------------------------------------------------
	// Create an FLTK window and display
	//------------------------------------------------------
	typedef Application<ImageType> ApplicationType;
	ApplicationType *mainWindow = new ApplicationType("Image Registration Tool");
	mainWindow->SetImageInputs(fixedImage, movingImage);
	mainWindow->show();
	Fl::run();

	CombineImages<ImageType, 
				  WriterType, 
				  AbsDiffFilterType>(fixedImage,
									  movingImage,
									  mainWindow->GetResultTransform(),
									  resultWriter);

	try
	{
		resultWriter->Update();
	}
  	catch (itk::ExceptionObject & e)
	{
    	std::cout << e << std::endl;
	    return EXIT_FAILURE;
	}

	return 0;
}
