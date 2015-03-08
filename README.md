# ITK_Haralick

ITK_Haralick is a plugin for the [ITK](http://www.itk.org/) Segmentation and Registration Toolkit that propose to compute moving Haralick texture features.

It is designed to be quick, and is thus independant from this [filter](http://www.itk.org/Doxygen/html/classitk_1_1Statistics_1_1HistogramToTextureFeaturesFilter.html) available in ITK. To ensure high performances, this filter only works on scalar images with a non-floating pixel type.

## How to use

This piece of code will compute moving Haralick texture features on a 3D unsigned char image using 16 grey levels, with a 5x5x5 window and for the 1,0,0 offset.

    #include <itkImage.h>
    #include <itkVectorImage.h>
    #include "itkScalarImageToHaralickTextureFeaturesImageFilter.h"

    typedef itk::Image< unsigned char, 3 > InputImageType;
    typedef itk::VectorImage< float, 3 > OutputImageType;

    typedef typename itk::Statistics::ScalarImageToHaralickTextureFeaturesImageFilter< InputImageType, typename OutputImageType::PixelType::ValueType > HaralickFilter;

    typename HaralickFilter::Pointer haralickImageComputer = HaralickFilter::New();
    haralickImageComputer->SetInput(your_rescaled_input_image);
    haralickImageComputer->SetNumberOfBinsPerAxis(16);

    typename HaralickFilter::RadiusType window_radius = {{5, 5, 5}};
    haralickImageComputer->SetWindowRadius(window_radius);

    typename HaralickFilter::OffsetVectorType::Pointer offsetV =
        HaralickFilter::OffsetVectorType::New();
    typename HaralickFilter::OffsetType offset;
    offset[0] = 1;
    offset[1] = 0;
    offset[2] = 0;

    offsetV->push_back(offset);

    haralickImageComputer->SetOffsets(offsetV);

    haralickImageComputer->Update();

    // haralickImageComputer->GetOutput();

The number of grey levels in your input image need to be reduce beforehand. This code will do it:

    #include <itkRescaleIntensityImageFilter.h>

    typedef itk::RescaleIntensityImageFilter< InputImageType, InputImageType > RescaleFilter;
    typename RescaleFilter::Pointer rescaler = RescaleFilter::New();
    rescaler->SetInput(input_image);
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(16 - 1);
    rescaler->Update();
    // Then use rescaler->GetOutput() as the input of the ScalarImageToHaralickTextureFeaturesImageFilter filter.

The `ScalarImageToHaralickTextureFeaturesImageFilter` filter is templated over the type of input image (the scalar image) and the type of the measurement that will be used to represent the Haralick features in the output image (a floating type, `float` is probably enough).

The output of the `ScalarImageToHaralickTextureFeaturesImageFilter` filter is an `itk::VectorImage` where each pixel is an `itk::VariableLengthVector` containing the following 10 Haralick normalized features: angular second moment, contrast, variance, inverse difference moment, sum average, sum variance, sum entropy, entropy, difference variance and difference entropy (other Haralick features have been left out due to numerical instability). Each feature is normalized from its native range to the [0; 1] range.

The `ScalarImageToHaralickTextureFeaturesImageFilter` filter allow you to use several offsets at the same time (the underlying cooccurrence matrix will be computed using all offsets), but it has never be really tested and the benefit of doing that has never been really proven.

## How to build

Run CMake, and during the "configure" step, set the `ITK_DIR` variable to the `lib/cmake/ITK-X.Y/` directory of your ITK installation. Also set the `CMAKE_INSTALL_DIR` variable.

Then, build and install the project.

## How to use it in your project

Add the next two lines to your `CMakeLists.txt`:

    find_package(ITK_Haralick REQUIRED)
    include_directories(${ITK_Haralick_INCLUDE_DIRS})

During the "configure" step, set the `ITK_Haralick_DIR` to the `lib/cmake/` directory where ITK_Haralick is installed.

No need to link anything, ITK_Haralick is an header only library.

Note: at the time it was written, it wasn't planned to make it an [ITK module](http://insightsoftwareconsortium.github.io/ITKBarCamp-doc/ITK/ConstructITKModule/index.html). If someone want to make it a module, your contribution is welcome.

## License

This tool is released under the terms of the MIT License. See the LICENSE.txt file for more details.