#define HAS_VTK 1
#define _IS_DEBUG 1

#include "vtkPointData.h"
#include <vtkPointPicker.h>
#include <vtkCommand.h>
#include <vtkMarchingCubes.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkCamera.h>
#include <vtkMarchingCubes.h>
#include <vtkVectorNorm.h>
#include <vtkDataSetMapper.h>
#include <vtkImageToPolyDataFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkLookupTable.h>
#include <vtkSphereSource.h>
#include <vtkCallbackCommand.h>
#include <vtkProperty.h>
#include <vtkImagePlaneWidget.h>
#include <vtkImageActor.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkPolyDataReader.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMaskPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointLocator.h>


#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;


void Get_Mean(char* source_poly_fn, char* target_poly_fn, char* output_poly_fn) 
{
    
    double target_scalar, source_scalar, mean_scalar;

  vtkSmartPointer<vtkPolyData> source_poly =vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> target_poly =vtkSmartPointer<vtkPolyData>::New();  
  
  // to search for closest point on polyWithColors for transfering the colors from that shell 
    
  vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New(); 
  reader1->SetFileName(source_poly_fn); 
  reader1->Update();
  source_poly = reader1->GetOutput();

  vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New(); 
  reader2->SetFileName(target_poly_fn); 
  reader2->Update();
  target_poly = reader2->GetOutput();


 vtkSmartPointer<vtkFloatArray> mean_scalars = vtkSmartPointer<vtkFloatArray>::New();
  for (vtkIdType i = 0; i < target_poly->GetNumberOfPoints(); ++i) {
       mean_scalars->InsertNextTuple1(0);  
  }

   vtkSmartPointer<vtkFloatArray> target_scalars = vtkSmartPointer<vtkFloatArray>::New();
    target_scalars = vtkFloatArray::SafeDownCast(target_poly->GetPointData()->GetScalars());

   vtkSmartPointer<vtkFloatArray> source_scalars = vtkSmartPointer<vtkFloatArray>::New();
    source_scalars = vtkFloatArray::SafeDownCast(source_poly->GetPointData()->GetScalars());
    

    // iterating through each point in source to find closest target point hit  
    for (vtkIdType i = 0; i < source_poly->GetNumberOfPoints(); ++i) {

        mean_scalar = 0;
        source_scalar = source_scalars->GetTuple1(i);
        target_scalar = target_scalars->GetTuple1(i);

        if (source_scalar > 0 && target_scalar > 0)
        {
            mean_scalar = (source_scalar+target_scalar) / 2; 
        }
        else if (source_scalar ==0)         // we assume 0 to be no data 
        {
            mean_scalar = target_scalar; 
        }
        else if (target_scalar == 0)       // we assume 0 to be no data
        {
            mean_scalar = source_scalar;
        }
        
        mean_scalars->SetTuple1(i, mean_scalar);

    }

    target_poly->GetPointData()->SetScalars(mean_scalars);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(target_poly);
    writer->SetFileName(output_poly_fn);
    writer->Update();

}



int main(int argc, char **argv)
{
  
    char* target_poly_fn, *source_poly_fn, *output_poly_fn, *output_txt_fn;

    vector<string> atlas_fns;
   
    if (argc != 4)
    {
        cerr << "Not enough parameters\nUsage: <shell_1.vtk> <shell_2.vtk> <mean.vtk>"<< endl; 
        exit(1);
    }

     source_poly_fn = argv[1];
    target_poly_fn = argv[2];
    output_poly_fn = argv[3];
    
    
    Get_Mean(source_poly_fn, target_poly_fn, output_poly_fn);
    

    
    
  
}