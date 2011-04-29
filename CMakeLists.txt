PROJECT(CP821)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(ITK REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
FIND_PACKAGE(FLTK REQUIRED)
FIND_PACKAGE(X11)

INCLUDE(${ITK_USE_FILE})
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(ImageReg itkImageToVTKImageFilter.txx
			itkRegistrationPipelineManager.cxx
			vtkFlRenderWindowInteractor.cxx
			vtkItkSourceViewer.cxx
			vtkPointSelectorStyle.cxx
			IterativeClosestPoint.cxx
		    PerformRegistrationButton.h
		    InputImageFiltersDialog.h)

INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIRS})

IF(X11_FOUND)
	TARGET_LINK_LIBRARIES(ImageReg ${X11_LIBRARIES} Xft Xinerama)
ENDIF(X11_FOUND)

SET(CMAKE_BUILD_TYPE "debug")

TARGET_LINK_LIBRARIES(ImageReg ITKCommon ITKStatistics ITKIO ITKNumerics 
                               vtkRendering vtkCommon vtkImaging vtkInfovis
			       fltk fltk_forms fltk_gl ITKBasicFilters)
