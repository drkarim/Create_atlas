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
#include <vtkCallbackCommand.h>u
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
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointLocator.h>


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm> 

using namespace std;


bool ReadFileNames(const char* fn, vector<string>& filename_list)
{
	ifstream infile(fn); 
	cout << "Reading filename list .. \n"; 
	std::string line;
	while (getline(infile, line))
	{
		istringstream iss(line);
		string filename; 
		if (!(iss >> filename)) { 
			cout << "Error reading file containing filenames, check format" << endl; 
			return false; 
		} // error
		else {
			cout << "Found " << filename << endl; 
			filename_list.push_back(filename);
		}

										 // process pair (a,b)
	}
	return true; 
}

void PrepareScalarContainer(string first_file, vector<vector< double>>& scalar_container)
{
	cout << "\nPreparing container structure with first file specified in list " << first_file << endl;
	vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
	reader1->SetFileName(first_file.c_str());
	reader1->Update();

	vtkSmartPointer<vtkPolyData> shell_poly = vtkSmartPointer<vtkPolyData>::New();
	shell_poly = reader1->GetOutput();

	vtkSmartPointer<vtkFloatArray> scalar_array = vtkSmartPointer<vtkFloatArray>::New();

	scalar_array = vtkFloatArray::SafeDownCast(shell_poly->GetPointData()->GetScalars());

	// construct the structure 
	/*
	for (int i = 0; i < shell_poly->GetNumberOfPoints(); ++i) {
		scalar_container.push_back(i);
	}*/
	
	// insert 
	cout << "Constructing ... "; 
	for (int i = 0; i < shell_poly->GetNumberOfPoints(); ++i) {
		
		double this_scalar = scalar_array->GetTuple1(i);
		vector<double> temp; 
		temp.push_back(this_scalar); 
		scalar_container.push_back(temp);
	}
	cout << "Finished prearing container structures!\n";

}

void InsertIntoContainer(string poly_data_fn, vector<vector< double>>& scalar_container)
{
	vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
	reader1->SetFileName(poly_data_fn.c_str());
	reader1->Update();

	cout << "Processing next file for atlas: " << poly_data_fn << " ... "; 
	vtkSmartPointer<vtkPolyData> shell_poly = vtkSmartPointer<vtkPolyData>::New();
	shell_poly = reader1->GetOutput();

	vtkSmartPointer<vtkFloatArray> scalar_array = vtkSmartPointer<vtkFloatArray>::New();
	scalar_array = vtkFloatArray::SafeDownCast(shell_poly->GetPointData()->GetScalars());

	for (vtkIdType i = 0; i < shell_poly->GetNumberOfPoints(); ++i) {

		double this_scalar = scalar_array->GetTuple1(i);
		scalar_container[i].push_back(this_scalar);
	}
	cout << "completed!" << endl;
}

// See Answer to https://stackoverflow.com/questions/2114797/compute-median-of-values-stored-in-vector-c

double CalcMedian(vector<double> scores)
{
	double median;
	size_t size = scores.size();

	sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}

double CalcMean(vector<double> scores)
{
	double sum = 0; double n = 0; 
	for (int i = 0; i < scores.size(); i++)
	{
		sum += scores[i];
		n++;
	}

	return sum / n; 
}

void CreateFinalAtlas(string first_file, const char* output_file, vector<vector< double>>& scalar_container, int mean_or_median)
{
	vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
	reader1->SetFileName(first_file.c_str());
	reader1->Update();

	vtkSmartPointer<vtkPolyData> shell_poly = vtkSmartPointer<vtkPolyData>::New();
	shell_poly = reader1->GetOutput();

	vtkSmartPointer<vtkFloatArray> atlas_scalars = vtkSmartPointer<vtkFloatArray>::New();
	atlas_scalars->SetNumberOfComponents(1); 

	for (vtkIdType i = 0; i < scalar_container.size(); ++i) {
		
		double atlas_value = -1; 
		if (mean_or_median == 1)
		{
			// computing mean 
			double mean = CalcMean(scalar_container[i]);
			atlas_value = mean; 
		}
		else if (mean_or_median == 2)
		{
			// compute median 
			double median = CalcMedian(scalar_container[i]); 
			atlas_value = median; 
		}
		

		atlas_scalars->InsertNextTuple1(atlas_value);
	}

	shell_poly->GetPointData()->SetScalars(atlas_scalars);

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(shell_poly);
	writer->SetFileName(output_file);
	writer->Update();

}


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
	bool foundArgs1 = false, foundArgs2 = false;
	int method = 1; // default is mean 
	const char* input_f, *output_f; 
	vector<string> filename_list;

    vector<string> atlas_fns;

	vector<vector<double> > scalar_container; 
	
	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++) {
			if (i + 1 != argc) {
				if (string(argv[i]) == "-i") {
					input_f = argv[i + 1];
					foundArgs1 = true;
				}
				else if (string(argv[i]) == "-m") {
					method = atoi(argv[i + 1]);
				}
				else if (string(argv[i]) == "-o") {
					output_f = argv[i + 1];
					foundArgs2 = true;
				}

			}
		}
	}
   
	if (!(foundArgs1 && foundArgs2))
	{
		cerr << "Cheeck your parameters\n\nUsage: \n\t-i <file_with_filenmes.txt> \n\t-m <method: 1=mean, 2=median> \n\t-o <atlas.vtk>" << endl;
		exit(1);
	}


   /*  source_poly_fn = argv[1];
    target_poly_fn = argv[2];
    output_poly_fn = argv[3];
    
    
    Get_Mean(source_poly_fn, target_poly_fn, output_poly_fn);*/ 

	if (!ReadFileNames(input_f, filename_list))
	{
		cerr << "Check the format of the file list. It should be \n\n\t\tfilename_1.vtk\n\t\tfilename.vtk\n\t\t... and so on" << endl;
		exit(1);
	}
	else {
		
		if (filename_list.size() > 0)
		{
			PrepareScalarContainer(filename_list[0], scalar_container);

			for (int i = 1; i < filename_list.size(); i++)
			{
				InsertIntoContainer(filename_list[i], scalar_container); 
			}

			CreateFinalAtlas(filename_list[0], output_f, scalar_container, method);

		}
		else {
			cerr << "File list empty, no VTK files to read for creating atlas" << endl; 
			exit(1); 
		}
	}
    

    
    
  
}