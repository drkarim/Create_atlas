# Creating an atlas from VTK files 
This repository contains code that will create an atlas from a list of [VTK](https://www.vtk.org/) files 

## Dependencies
Note that this program uses [VTK - Visualisation Toolkit C++ Library](https://www.vtk.org/) 

## Usage 
The usage for the tool ```create_atlas``` is as follows:
```
create_atlas
     -i <list of filenames with extension in each new line, as a text file> 
     -m <1=mean or 2=median for creating atlas> 
     -o <output atlas file name>
```

## Important notes
it uses the first file listed in the filename list (```-i``` switch) as the structure of the atlas. The structure here refers to both the spatial location and number of vertices. 

The output of the atlas is same as the input and is in VTK format. 

## Author 
```
Dr. Rashed Karim 
Department of Biomedical Engineering 
King's College London 
```
