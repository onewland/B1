#ifndef _itk_registration_pipeline_manager_cxx
#define _itk_registration_pipeline_manager_cxx
#include "itkRegistrationPipelineManager.h"

namespace itk
{

template<class TInputImage, int Dimension>
RegistrationPipelineManager<TInputImage,Dimension>::
RegistrationPipelineManager() :
	numberOfIterations(500),
	gradientTolerance(1E-5),
	valueTolerance(1E-5),
	epsilonFunction(1E-6)
{
	transform = TransformType::New();
	optimizer = OptimizerType::New();
	optimizer->SetUseCostFunctionGradient(false);
	registration = RegistrationType::New();
	metric = MetricType::New();
}

template<class TInputImage, int Dimension>
void
RegistrationPipelineManager<TInputImage,Dimension>::
SetRegistrationPoints(PointSetTypePointer fixedPoints, 
        			   PointSetTypePointer movingPoints)
{
	fixedPointSet = fixedPoints;
	movingPointSet = movingPoints;
}

template<class TInputImage, int Dimension>
void
RegistrationPipelineManager<TInputImage,Dimension>::
PerformRegistration()
{
	double guess_angle = 3*Math::pi/2;
	OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
	scales[0] = Math::pi/1000.0;
	scales[1] = 1.0/10.0;
	scales[2] = 1.0/10.0;
	scales[3] = 1.0/100.0;
	scales[4] = 1.0/100.0;
	
	ParametersType defaults(transform->GetNumberOfParameters());

	// rotation
	defaults.SetElement(0, guess_angle);
	// center x,y
	defaults.SetElement(1, 0.0);
	defaults.SetElement(2, 0.0);
	// translation x,y 
	defaults.SetElement(3, 0.0);
	defaults.SetElement(4, 0.0);
	transform->SetParameters(defaults);

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

	registration->StartRegistration();
}

template<class TInputImage, int Dimension>
itk::CenteredRigid2DTransform<double>::Pointer
RegistrationPipelineManager<TInputImage,Dimension>::
GetTransform()
{
	return transform;
}


}
#endif
