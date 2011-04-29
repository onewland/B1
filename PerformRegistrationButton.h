#ifndef __perform_registration_button_h
#define __perform_registration_button_h

#include <FL/Fl_Button.H>
#include <vtkArrayIterator.h>
#include "vtkItkSourceViewer.h"
#include "itkRegistrationPipelineManager.h"

typedef itk::CenteredRigid2DTransform< double > TransformType;

template<class TInputImage>
class PerformRegistrationButton : public Fl_Button {
private:
	static const unsigned int Dimension = 2;
	typedef itk::PointSet < float, Dimension > PointSetType;
	typedef TInputImage ImageType;
	typedef vtkITKSourceViewer<ImageType> InputViewer;
	typedef PointSetType::PointType PointType;
	typedef PointSetType::PointsContainer PointsContainer;
	typedef PointsContainer::Iterator PointsIterator;
	typedef itk::RegistrationPipelineManager<ImageType, Dimension> RegistrationManager;
	typedef typename RegistrationManager::Pointer RegistrationManagerPointer;
	typedef PerformRegistrationButton<ImageType> Self;
	TransformType::Pointer transform;

	InputViewer *fixedPointSource;
	InputViewer *movingPointSource;

public:
	PerformRegistrationButton(int x, int y, int w, int h, 
		const char *s = "Perform Registration") : Fl_Button(x,y,w,h,s)
	{
		type(FL_NORMAL_BUTTON);
		callback(PerformRegistrationCallback);
	}

	virtual TransformType::Pointer GetTransform()
	{
		return transform;
	}

	virtual void SetFixedPointSource(InputViewer *fps)
	{
		fixedPointSource = fps;
	}

	virtual void SetMovingPointSource(InputViewer *mps)
	{
		movingPointSource = mps;
	}

	static void PerformRegistrationCallback(Fl_Widget *registration_button)
	{
		Self *pbr = reinterpret_cast<Self *>(registration_button);
		pbr->Perform();
	}

	virtual void Perform() 
	{
		vtkSmartPointer<vtkPoints> fixedVTKPoints = fixedPointSource->GetSelectedPoints();
		vtkSmartPointer<vtkPoints> movingVTKPoints = movingPointSource->GetSelectedPoints();
		PointsContainer::Pointer fixedPointContainer = PointsContainer::New();
		PointsContainer::Pointer movingPointContainer = PointsContainer::New();
		PointSetType::Pointer fixedPointSet = PointSetType::New();
		PointSetType::Pointer movingPointSet = PointSetType::New();
		RegistrationManagerPointer registrationManager = RegistrationManager::New();
		unsigned int pointId = 0;
		unsigned int fixedCount = fixedVTKPoints->GetNumberOfPoints();
		unsigned int movingCount = movingVTKPoints->GetNumberOfPoints();
		for(pointId = 0; pointId < fixedCount; pointId++)
		{
			double* fixedPoint = fixedVTKPoints->GetPoint(pointId);
			float fixedPointF[] = { fixedPoint[0], 500 - fixedPoint[1] };
			std::cout << "Fixed point " << pointId << " = ( " 
					  << fixedPointF[0] << ", " 
					  << fixedPointF[1] << ")" << std::endl;
			PointType p(fixedPointF);
			fixedPointContainer->InsertElement(pointId, p);
		}
		for(pointId = 0; pointId < movingCount; pointId++)
		{
			double* movingPoint = movingVTKPoints->GetPoint(pointId);
			float movingPointF[] = { movingPoint[0], 500 - movingPoint[1] };
			std::cout << "Moving point " << pointId << " = ( " 
					  << movingPointF[0] << ", " 
					  << movingPointF[1] << ")" << std::endl;

			PointType p(movingPointF);
			movingPointContainer->InsertElement(pointId, p);
		}

		fixedPointContainer->Print(std::cout);
		movingPointContainer->Print(std::cout);
		fixedPointSet->SetPoints(fixedPointContainer);
		movingPointSet->SetPoints(movingPointContainer);
		registrationManager->SetRegistrationPoints(fixedPointSet, movingPointSet);
		registrationManager->PerformRegistration();
		std::cout << "Solution = " 
				  << registrationManager->GetTransform()->GetParameters() 
				  << std::endl;
		transform = registrationManager->GetTransform();
	}
};

#endif
