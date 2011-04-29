#ifndef _itk_registration_pipeline_manager_h
#define _itk_registration_pipeline_manager_h

#include <itkLevenbergMarquardtOptimizer.h>
#include <itkCommand.h>
#include <itkPointSet.h>
#include <itkEuclideanDistancePointMetric.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkPointSetToPointSetRegistrationMethod.h>

// Use this class to debug registration
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

namespace itk 
{

template<class TInputImage, int Dimension>
class RegistrationPipelineManager : public itk::Object
{
public:
	typedef RegistrationPipelineManager Self;
	typedef SmartPointer<Self> Pointer;
	typedef SmartPointer<const Self> ConstPointer;

	typedef itk::PointSet< float, Dimension > PointSetType;
	typedef typename PointSetType::Pointer PointSetTypePointer;
	typedef typename PointSetType::PointsContainer PointsContainer;
	typedef TInputImage InputImageType;
	typedef typename itk::EuclideanDistancePointMetric<PointSetType, PointSetType>
	    MetricType;
	typedef typename MetricType::Pointer MetricTypePointer; 
	typedef typename MetricType::TransformType TransformBaseType;
  	typedef typename TransformBaseType::ParametersType ParametersType;
	typedef typename TransformBaseType::JacobianType JacobianType;

	typedef typename itk::CenteredRigid2DTransform< double > TransformType;
	typedef typename itk::LevenbergMarquardtOptimizer OptimizerType;
	typedef typename itk::PointSetToPointSetRegistrationMethod
		<PointSetType, PointSetType> RegistrationType;
	typedef typename RegistrationType::Pointer RegistrationTypePointer;

	itkNewMacro(Self);

	void SetRegistrationPoints(PointSetTypePointer fixedPoints, 
						       PointSetTypePointer movingPoints);
	void PerformRegistration();
	TransformType::Pointer GetTransform();

protected:
	RegistrationPipelineManager();

private:
	PointSetTypePointer fixedPointSet;// = PointSetType::New();
	PointSetTypePointer movingPointSet;// = PointSetType::New();
	TransformType::Pointer transform;
	OptimizerType::Pointer optimizer;
	RegistrationTypePointer registration;
	MetricTypePointer metric;

	const unsigned long numberOfIterations;// = 500;
	const double gradientTolerance;// = 1e-5;	// convergence criterion
	const double valueTolerance;// = 1e-5;	// convergence criterion
	const double epsilonFunction;// = 1e-6;	// convergence criterion
};

}

#include "itkRegistrationPipelineManager.cxx"
#endif
