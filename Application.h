#ifndef __application_h
#define __application_h

#include "vtkItkSourceViewer.h"
#include "PerformRegistrationButton.h"

template<class TInputImage>
class Application : public Fl_Window {
private:
	typedef TInputImage ImageType;
	typedef vtkITKSourceViewer<ImageType> InputViewer;
	typedef PerformRegistrationButton<ImageType> RegisterButton;
	typedef Application<ImageType> Self;
	InputViewer *fl_vtk_viewer_left;
	InputViewer *fl_vtk_viewer_right;
	RegisterButton *register_button;
	Fl_Button *show_preprocessor_fixed;
	Fl_Button *show_preprocessor_moving;

	static void ShowFixedOptionsCallback(Fl_Widget *widget, void *self)
	{
		Self *app = reinterpret_cast<Self *>(self);
		app->ShowFixedPreprocessorWindow();
	}

	static void ShowMovingOptionsCallback(Fl_Widget *widget, void *self)
	{
		Self *app = reinterpret_cast<Self *>(self);
		app->ShowMovingPreprocessorWindow();
	}


public:
	Application(const char *windowTitle) :
		Fl_Window(1010, 600, windowTitle)
	{
		fl_vtk_viewer_left = new InputViewer(0,0,500,500);
		fl_vtk_viewer_right = new InputViewer(510,0,500,500);
		fl_vtk_viewer_left->Initialize();
		fl_vtk_viewer_right->Initialize();
		register_button = new RegisterButton(25,510,250,50);
		register_button->SetFixedPointSource(fl_vtk_viewer_left);
		register_button->SetMovingPointSource(fl_vtk_viewer_right);
		show_preprocessor_fixed = new Fl_Button(300,510,200,50, "Fixed Image Options");
		show_preprocessor_fixed->callback(ShowFixedOptionsCallback, this);
		show_preprocessor_moving = new Fl_Button(525,510,200,50, "Moving Image Options");
		show_preprocessor_moving->callback(ShowMovingOptionsCallback, this);
		end();
		
		// these must be done after end so they don't get added
		// to the main window
		fl_vtk_viewer_left->MakeFilterDialog("Fixed Pre-Processing Options");
		fl_vtk_viewer_right->MakeFilterDialog("Moving Pre-Processing Options");
	}

	void SetImageInputs(const ImageType *fixedImage, 
						const ImageType *movingImage)
	{
		fl_vtk_viewer_left->SetImageInput(fixedImage);
		fl_vtk_viewer_right->SetImageInput(movingImage);
	}
	
	virtual void ShowFixedPreprocessorWindow()
	{
		fl_vtk_viewer_left->GetFilterDialog()->show();
	}

	virtual void ShowMovingPreprocessorWindow()
	{
		fl_vtk_viewer_right->GetFilterDialog()->show();
	}

	virtual void show()
	{
		Fl_Window::show();
		fl_vtk_viewer_left->show();
		fl_vtk_viewer_right->show();
	}

	virtual TransformType::Pointer GetResultTransform()
	{
		return register_button->GetTransform();
	}
};

#endif
