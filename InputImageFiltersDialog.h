#include <FL/Fl_Window.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include "vtkItkSourceViewer.h"

template <class TInputImage> class vtkITKSourceViewer;

template <class TInputImage> 
class InputImageFiltersDialog : public Fl_Window {
	typedef TInputImage InputImageType;
	typedef InputImageFiltersDialog<InputImageType> Self;
	typedef vtkITKSourceViewer<TInputImage> InputViewer;

	Fl_Slider *intensityWindowMinSlider;
	Fl_Slider *intensityWindowMaxSlider;
	Fl_Slider *intensityOffsetSlider;
	Fl_Group *intensityWindowBox;
	Fl_Button *intensityWindowResetButton;
	InputViewer *imageInputViewer;

	static void MinSliderCallback(Fl_Widget *widget, void *data)
	{
		Fl_Slider *slider = reinterpret_cast<Fl_Slider*>(widget);
		Self *dialog = reinterpret_cast<Self *>(data);
		dialog->UpdateWindowIntensityMin(slider->value());
	}

	static void MaxSliderCallback(Fl_Widget *widget, void *data)
	{
		Fl_Slider *slider = reinterpret_cast<Fl_Slider*>(widget);
		Self *dialog = reinterpret_cast<Self *>(data);
		dialog->UpdateWindowIntensityMax(slider->value());
	}

	static void ResetIntensityWindowCallback(Fl_Widget *widget, void *data)
	{
		Self *dialog = reinterpret_cast<Self *>(data);
		dialog->intensityWindowMinSlider->value(0);
		dialog->intensityWindowMaxSlider->value(255);
		dialog->UpdateWindowIntensityMin(0);
		dialog->UpdateWindowIntensityMax(255);
	}

	static void OffsetSliderCallback(Fl_Widget *widget, void *data)
	{
		Fl_Slider *slider = reinterpret_cast<Fl_Slider*>(widget);
		Self *dialog = reinterpret_cast<Self *>(data);
		dialog->UpdateIntensityOffset(slider->value());
	}

	void UpdateWindowIntensityMin(double intensity)
	{
		imageInputViewer->SetIntensityWindowMin(intensity);
	}

	void UpdateWindowIntensityMax(double intensity)
	{
		imageInputViewer->SetIntensityWindowMax(intensity);
	}

	void UpdateIntensityOffset(double intensity)
	{
		imageInputViewer->SetIntensityOffset(intensity);
	}

	void CreateIntensityWindowBox() 
	{
		intensityWindowBox = new Fl_Group(10, 20, 256, 120, "Intensity Window");
		intensityWindowBox->align(FL_ALIGN_TOP);
		intensityWindowBox->box(FL_BORDER_BOX);
		{
			int leftOffset = 15;
			intensityWindowMinSlider = 
				new Fl_Slider(leftOffset,30,100,30,"Min Intensity");
			intensityWindowMinSlider->type(FL_HOR_NICE_SLIDER);
			intensityWindowMinSlider->range(0,254);
			intensityWindowMinSlider->value(0);
			intensityWindowMinSlider->align(FL_ALIGN_RIGHT);
			intensityWindowMinSlider->callback(MinSliderCallback,this);
			intensityWindowMaxSlider = 
				new Fl_Slider(leftOffset,60,100,30,"Max Intensity");
			intensityWindowMaxSlider->type(FL_HOR_NICE_SLIDER);
			intensityWindowMaxSlider->range(1,255);
			intensityWindowMaxSlider->value(255);
			intensityWindowMaxSlider->align(FL_ALIGN_RIGHT);
			intensityWindowMaxSlider->callback(MaxSliderCallback,this);
			intensityWindowResetButton = 
				new Fl_Button(leftOffset,100,244,20,"Reset Intensity Window");
			intensityWindowResetButton->callback(ResetIntensityWindowCallback,this);
		}
		intensityWindowBox->end();
	}
public:
	InputImageFiltersDialog(const char *title) 
		: Fl_Window(500,300, title) 
	{
		CreateIntensityWindowBox();

		intensityOffsetSlider = new Fl_Slider(10, 160, 256, 30, 
			"Intensity Offset (Brightness)");
		intensityOffsetSlider->type(FL_HOR_NICE_SLIDER);
		intensityOffsetSlider->range(0,255);
		intensityOffsetSlider->value(0);
		intensityOffsetSlider->align(FL_ALIGN_RIGHT);
		intensityOffsetSlider->callback(OffsetSliderCallback, this);
		
		end();
	}

	void SetInputViewer(InputViewer *viewer)
	{
		imageInputViewer = viewer;
	}
};
