import os
import subprocess
import shutil


#%%%%%%%%%%%%%%%%% Start of Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunDir = 'RAE2822'
GridDir= 'CreateBlocks/grid'
NumberOfBlocks = 6

def SetInput(Control, Scheme, Flow, OutputControl, ResidualControl):
    Control['CFL'] = 100.0
    Control['LoadLevel'] = 0
    Control['MaxIterations'] = 10000
    Control['SaveIterations'] = 100
    Control['OutputFileFormat'] = 'tecplot'
    Control['OutputDataFormat'] = 'ASCII'
    Control['InputFileFormat'] = 'tecplot'
    Control['InputDataFormat'] = 'ASCII'
    Control['Precision'] = 6
    Control['Purge'] = 2
    Control['ResidualWriteInterval'] = 100
    Control['Tolerance'] = "1e-13 Continuity_abs"
    Control['DebugLevel'] = 5
    
    Scheme['InviscidFlux'] = 'ausm'
    Scheme['FaceState'] = 'muscl'
    Scheme['Limiter'] = '0 0 0  0 0 0'
    Scheme['TurbulenceLimiter'] = '1 1 1'
    Scheme['TurbulenceModel']='sst'
    Scheme['TransitionModel']='none'
    Scheme['TimeStep']='l'
    Scheme['TimeIntegration']='plusgs'
    Scheme['HigherOrderBC']='0'
    
    Flow["NumberOfVariables"] = 5
    Flow["DensityInf"] = 1.2
    Flow["UInf"] = 252.9
    Flow["VInf"] = 10.2
    Flow["WInf"] = 0.0
    Flow["PressureInf"] = 103338.0
    Flow["TurbulenceIntensity"] = 0.03873
    Flow["ViscosityRatio"] = 0.01
    Flow["Intermittency"] = 1.0
    Flow["ReferenceViscosity"] = 1.4243e-5
    Flow["ViscosityLaw"] = "sutherland_law"
    Flow["ReferenceTemp"] = 300
    Flow["SutherlandTemp"] = 110.5
    Flow["PrandtlNumbers"] = "0.72 0.9"
    Flow["SpecificHeatRatio"]=1.4
    Flow["GasConstant"]=287.0

    OutputControl['Out'] = ["Velocity", "Density", "Pressure", "Mu"]
    OutputControl['In'] = ["Velocity", "Density", "Pressure", "Mu"]
    ResidualControl['Out'] = ["Mass_abs", "Viscous_abs", "Continuity_abs"]
    BoundaryConditions = [-4, -4, -5, -8, -6, -6]
    return BoundaryConditions

#%%%%%%%%%%%%%%%%%%%% End of Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BC={-1:'SUPERSONIC INFLOW (DIRICHLET)', -2:'SUPERSONIC OUTFLOW (EXTRAPOLATION)', -3:'SUBSONIC INFLOW (MASS-FLOW RATE FIXED)', -4:'SUBSONIC OUTFLOW (PRESSURE FIXED)', -5:'WALL (NO SLIP)', -6:'SYMMETRY', -7:'POLE', -8:'FAR-FIELD', -11:'TOTAL INLET'}


def SetExpectedInput(ExpectedControl, ExpectedScheme, ExpectedFlow, ExpectedOutputControl, ExpectedResidualControl):
    ExpectedControl['CFL'] = ['0.25', '0.5', '0.75', '1.0', '10', '100']
    ExpectedControl['LoadLevel'] = '10'
    ExpectedControl['MaxIterations'] = '1000000'
    ExpectedControl['SaveIterations'] = '100000'
    ExpectedControl['InputFileFormat'] = ['vtk', 'tecplot']
    ExpectedControl['InputDataFormat'] = ['ASCII']
    ExpectedControl['OutputFileFormat'] = ['vtk', 'tecplot']
    ExpectedControl['OutputDataFormat'] = ['ASCII']
    ExpectedControl['Precision'] = 6
    ExpectedControl['Purge'] = 2
    ExpectedControl['ResidualWriteInterval'] = 100
    ExpectedControl['Tolerance'] = "1e-13 Continuity_abs"
    ExpectedControl['DebugLevel'] = [1, 2, 3, 4, 5]
    ExpectedScheme['InviscidFlux'] = ['ausm', 'ldfss0', 'slau', 'van_leer', 'ausmUP', 'ausmP']
    ExpectedScheme['FaceState'] = ['muscl', 'none', 'ppm', 'weno']
    ExpectedScheme['Limiter'] = '1 1 1  0 0 0'
    ExpectedScheme['TurbulenceLimiter'] = '1 1 1'
    ExpectedScheme['TurbulenceModel']=['none', 'sst', 'kkl', 'sa', 'sst2003']
    ExpectedScheme['TransitionModel']=['none', 'bc', 'lctm2015']
    ExpectedScheme['TimeStep']=['g', 'l']
    ExpectedScheme['TimeIntegration']=['RK4', 'RK2', 'TVDRK2', 'TVDRK3', 'implicit', 'none', 'plusgs']
    ExpectedScheme['HigherOrderBC']='0'
    ExpectedFlow['ViscosityLaw'] = ['sutherland_law', 'constant']
    ExpectedOutputControl['Out']=['Velocity', 
                                 'Density', 
                                 'Pressure',
                                 'Mu',
                                 'Mu_t',
                                 'TKE',
                                 'Omega',
                                 'kL',
                                 'tv',
                                 'Wall_distance'
                                 'resnorm',
                                 'TKE_residue',
                                 'Mass_residue'
                                 'X_mom_residue',
                                 'Y_mom_residue',
                                 'Z_mom_residue',
                                 'energy_residue',
                                 'DuDx', 'Dudy', 'DuDz',
                                 'DvDx', 'DvDy', 'DvDz',
                                 'DwDx', 'DWDy', 'DwDz',
                                 'DTDx', 'DTDy', 'DTDz',
                                 'DtkDx', 'DtkDy', 'DtkDz',
                                 'DtwDx', 'DtwDy', 'DtwDz',
                                 'DtvDx', 'DtvDy', 'DtvDz',
                                 'DtkLDx', 'DtkLDy', 'DtkLDz',
                                 ]
    
    ExpectedOutputControl['In']=['Velocity', 
                                  'Density', 
                                  'Pressure',
                                  'viscosity',
                                  'TKE',
                                  'Omega',
                                  'kL',
                                  'tv',
                                  ]
    
    ExpectedResidualControl['Out'] = ["Mass_abs",
                                      "Viscous_abs",
                                      "Mass_abs",
                                      "Resnorm_abs",
                                      "Viscous_abs",
                                      "Turbulent_abs",
                                      "Continuity_abs",
                                      "X-mom_abs",
                                      "Y-mom_abs",
                                      "Z-mom_abs",
                                      "Energy_abs",
                                      "Mass_rel",
                                      "Resnorm_rel",
                                      "Viscous_rel",
                                      "Turublent_rel",
                                      "Continuity_abs",
                                      "X-mom_rel",
                                      "Y-mom_rel",
                                      "Z-mom_rel",
                                      "Energy_rel",
                                      "TKE_abs",
                                      "Tv_abs",
                                      "Dissipation_abs",
                                      "Omega_abs",
                                      "Kl_abs",
                                      "TKE_rel",
                                      "Tv_rel",
                                      "Dissipation_rel",
                                      "Omega_rel",
                                      "Kl_rel"
                                     ]

def CheckInput(ExpectedControl, ExpectedScheme, ExpectedFlow, ExpectedOutputControl, ExpectedResidualControl, Control, Scheme, Flow, OutputControl, ResidualControl):
    assert os.environ.get("FEST3D") is not None, 'FEST3D enviornment variable is not set\nPlease use following command to export the FEST3D variable and then rerun the python scirpt\n \n$export FEST3D="absolute/path/to/bin/FEST3D"\n'
    assert (Control['CFL'] > 0)
    assert (type(Control['LoadLevel']) == int and Control['LoadLevel'] >= 0)
    assert (type(Control['MaxIterations']) == int and Control['MaxIterations'] >= 0)
    assert (type(Control['SaveIterations']) == int and Control['SaveIterations'] >= 0)
    assert Control['SaveIterations'] <= Control['MaxIterations']
    assert Control['MaxIterations'] % Control['SaveIterations'] == 0 , "Change checkpoint interations; Likely to miss last Maxiteraion point."
    assert Control['InputFileFormat'] in ExpectedControl['InputFileFormat']
    assert Control['InputDataFormat'] in ExpectedControl['InputDataFormat']
    assert Control['OutputFileFormat'] in ExpectedControl['OutputFileFormat']
    assert Control['OutputDataFormat'] in ExpectedControl['OutputDataFormat']
    assert (type(Control['Precision']) == int and Control['Precision'] >= 4)
    assert (type(Control['Purge']) == int and Control['Purge'] >= 0)
    assert (type(Control['ResidualWriteInterval']) == int and Control['ResidualWriteInterval'] >=1 )
    assert Control['Tolerance'].split()[-1] in ResidualControl['Out']
    assert Control['DebugLevel'] in ExpectedControl['DebugLevel']
    assert Scheme['InviscidFlux'] in ExpectedScheme['InviscidFlux']
    assert Scheme['FaceState'] in ExpectedScheme['FaceState']
    assert Scheme['TurbulenceModel'] in ExpectedScheme['TurbulenceModel']
    assert Scheme['TransitionModel'] in ExpectedScheme['TransitionModel']
    assert Scheme['TimeIntegration'] in ExpectedScheme['TimeIntegration']
    assert Scheme['TimeStep'].split()[0] in ExpectedScheme['TimeStep']
    assert Scheme['HigherOrderBC'] in ExpectedScheme['HigherOrderBC']
    assert Flow['ViscosityLaw'] in ExpectedFlow['ViscosityLaw']
    assert all(variable in ExpectedOutputControl['Out'] for variable in OutputControl['Out'])
    assert all(variable in ExpectedOutputControl['Out'] for variable in OutputControl['In'])
    assert all(variable in ExpectedResidualControl['Out'] for variable in ResidualControl['Out'])
    #Number of grid files should be equal number of blocks as input
    assert len(next(os.walk(GridDir))[2]) == NumberOfBlocks, "Please run the Makefile in CreateBlocks directory to generate grids"
    
    

def WriteInputFiles(RootDir):
    WriteControlFile(Control, RootDir)
    WriteFvSchemeFile(Scheme, RootDir)
    WriteFlowFile(Flow, RootDir)
    WriteOutputControlFile(OutputControl, RootDir)
    WriteResidualControlFile(ResidualControl, RootDir)
    WriteStopFile('0', RootDir)



def WriteControlFile(Control, RootDir):
    with open(RootDir+"/system/control.md", "w+") as file:
        file.write("\n")
        file.write("Control File\n")
        file.write("===========\n")
        file.write("## CFL\n")
        file.write(str(Control['CFL'])+"\n\n")
        file.write("## Restart level\n")
        file.write(str(Control['LoadLevel'])+"\n\n")
        file.write("## Maximum Iterations\n")
        file.write(str(Control['MaxIterations'])+"\n\n")
        file.write("## Save Iterations\n")
        file.write(str(Control['SaveIterations'])+"\n\n")
        file.write("## Output File format\n")
        file.write(str(Control['OutputFileFormat'])+"\n\n")
        file.write("## Output Data format\n")
        file.write(str(Control['OutputDataFormat'])+"\n\n")
        file.write("## Input File format\n")
        file.write(str(Control['InputFileFormat'])+"\n\n")
        file.write("## Input Data format\n")
        file.write(str(Control['InputDataFormat'])+"\n\n")
        file.write("## Write Precision\n")
        file.write(str(Control['Precision'])+"\n\n")
        file.write("## Purge Write\n")
        file.write(str(Control['Purge'])+"\n\n")
        file.write("## Residual write Interval\n")
        file.write(str(Control['ResidualWriteInterval'])+"\n\n")
        file.write("## Tolerance\n")
        file.write(str(Control['Tolerance'])+"\n\n")
        file.write("## Debug Level\n")
        file.write(str(Control['DebugLevel'])+"\n\n")



def WriteFvSchemeFile(Scheme, RootDir):
    with open(RootDir+"/system/fvscheme.md", "w+") as file:
        file.write("\n")
        file.write("Scheme File\n")
        file.write("===========\n")
        file.write("## Inviscid Flux Scheme\n")
        file.write(str(Scheme['InviscidFlux'])+"\n\n")
        file.write("## Higher Order Method\n")
        file.write(str(Scheme['FaceState'])+"\n\n")
        file.write("## Switch: Limiter - Pressure based switch\n")
        file.write(str(Scheme['Limiter'])+"\n\n")
        file.write("## Turbulence Limiter Switch\n")
        file.write(str(Scheme['TurbulenceLimiter'])+"\n\n")
        file.write("## Turbulence model\n")
        file.write(str(Scheme['TurbulenceModel'])+"\n\n")
        file.write("## Transition model\n")
        file.write(str(Scheme['TransitionModel'])+"\n\n")
        file.write("## Time Step\n")
        file.write(str(Scheme['TimeStep'])+"\n\n")
        file.write("## Time Integration Method\n")
        file.write(str(Scheme['TimeIntegration'])+"\n\n")
        file.write("## Higher Order Boundary Conditions\n")
        file.write(str(Scheme['HigherOrderBC'])+"\n\n")



def WriteFlowFile(Flow, RootDir):
    with open(RootDir+"/system/flow.md", "w+") as file:
        file.write("\n")
        file.write("Flow File\n")
        file.write("===========\n")
        file.write("## Number of Variables\n")
        file.write(str(Flow['NumberOfVariables'])+"\n\n")
        file.write("## Free Stream Density\n")
        file.write(str(Flow['DensityInf'])+"\n\n")
        file.write("## Free Stream X-Speed\n")
        file.write(str(Flow['UInf'])+"\n\n")
        file.write("## Free Stream Y-Speed\n")
        file.write(str(Flow['VInf'])+"\n\n")
        file.write("## Free Stream Z-Speed\n")
        file.write(str(Flow['WInf'])+"\n\n")
        file.write("## Free Stream Pressure\n")
        file.write(str(Flow['PressureInf'])+"\n\n")
        file.write("## Free Stream Turbulence Intensity\n")
        file.write(str(Flow['TurbulenceIntensity'])+"\n\n")
        file.write("## Free Stream Viscosity Ratio\n")
        file.write(str(Flow['ViscosityRatio'])+"\n\n")
        file.write("## Free Stream Intermittency\n")
        file.write(str(Flow['Intermittency'])+"\n\n")
        file.write("## Reference Viscosity\n")
        file.write(str(Flow['ReferenceViscosity'])+"\n\n")
        file.write("## Viscosity law\n")
        file.write(str(Flow['ViscosityLaw'])+"\n\n")
        file.write("## Reference Temperature\n")
        file.write(str(Flow['ReferenceTemp'])+"\n\n")
        file.write("## Sutherland Temperature\n")
        file.write(str(Flow['SutherlandTemp'])+"\n\n")
        file.write("## Prandtl numbers\n")
        file.write(str(Flow['PrandtlNumbers'])+"\n\n")
        file.write("## Specific heat ratio\n")
        file.write(str(Flow['SpecificHeatRatio'])+"\n\n")
        file.write("## Gas Constant\n")
        file.write(str(Flow['GasConstant'])+"\n\n")



def WriteOutputControlFile(OutputControl, RootDir):
    with open(RootDir+"/system/output_control.md", "w+") as file:
        file.write("{\n")
        for item in  OutputControl['Out']:
            file.write("  "+item+"\n")
        file.write("}\n")
        file.write("{\n")
        for item in  OutputControl['In']:
            file.write("  "+item+"\n")
        file.write("}\n")


def WriteResidualControlFile(ResidualControl, RootDir):
    with open(RootDir+"/system/res_control.md", "w+") as file:
        file.write("{\n")
        for item in  ResidualControl['Out']:
            file.write("  "+item+"\n")
        file.write("}\n")

def WriteStopFile(Data, RootDir):
    with open(RootDir+"/system/stopfile", "w+") as file:
        file.write(Data)


def Create(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


ResidualControl = {}
OutputControl = {}
Scheme={}
Flow={}
Control={}
ExpectedScheme={}
ExpectedFlow={}
ExpectedControl={}
ExpectedResidualControl={}
ExpectedOutputControl={}
BoundaryConditions=[]

BoundaryConditions = SetInput(Control, Scheme, Flow, OutputControl, ResidualControl)
SetExpectedInput(ExpectedControl, ExpectedScheme, ExpectedFlow, ExpectedOutputControl, ExpectedResidualControl)
CheckInput(ExpectedControl, ExpectedScheme, ExpectedFlow, ExpectedOutputControl, ExpectedResidualControl, Control, Scheme, Flow, OutputControl, ResidualControl)

Create(RunDir)
Create(RunDir+'/system')
Create(RunDir+'/system/mesh')
Create(RunDir+'/system/mesh/gridfiles')
Create(RunDir+'/system/mesh/layout')
Create(RunDir+'/system/mesh/bc')
Create(RunDir+'/time_directories')
Create(RunDir+'/time_directories/aux')
Create(RunDir+'/pp')
Create(RunDir+'/pre')
Create(RunDir+'/bin')
WriteInputFiles(RunDir)

BCListFilename="/system/mesh/layout/bc_list"
BCValueFilename="/system/mesh/layout/bc_values"
GenBCFilename="/system/mesh/layout/generate_bc.cpp"
GenLayoutFilename="/system/mesh/layout/generate_layout.cpp"
CompilerFilename="/system/mesh/layout/compile.sh"
NameFillerFilename="/system/mesh/layout/fill_names.sh"
BashRunFilename="/run.sh"
GnuplotFilename="/gnplt"
VTKFillerFilename="/fill_vtk_name.sh"
GenBCFile="generate_bc.cpp"
GenLayoutFile="generate_layout.cpp"
CompileFile="compile.sh"
FillerFile="fill_names.sh"
BashRunFile="run.sh"
GnuplotFile="gnplt"
VTKFillerFile="fill_vtk_name.sh"

def WriteFullFile(RootDir,Filename, Data):
    with open(RootDir+Filename, "w+") as file:
        file.write(Data)

BCList = '''## BC LIST, BCs SEPERATED BY HYPHEN
## DON'T CHANGE INTERFACE VALUE
## SUPERSONIC INFLOW
-1 1-2-3-4-5
## SUPERSONIC OUTFLOW
-2 6-7-8-9-10
## SUBSONIC INFLOW
-3 1-2-3-4-10
## SUBSONIC OUTFLOW
-4 6-7-8-9-5
## WALL
-5 6-10-12-17
## SLIP PLANE
-6 6-10-11
## Pole
-7 15
## Far-field
-8 16
## Total inlet
-9 19
## INTERFACE
1 13
'''

BCValue='''## NUMBER BC_VALUE
 1 FIX_DENSITY
 2 FIX_X_SPEED
 3 FIX_Y_SPEED
 4 FIX_Z_SPEED
 5 FIX_PRESSURE
 6 COPY_DENSITY
 7 COPY_X_SPEED
 8 COPY_Y_SPEED
 9 COPY_Z_SPEED
 10 COPY_PRESSURE
 11 FLOW_TANGENCY
 12 NO_SLIP
 13 INTERFACE
 14 PERIODIC
 15 POLE
 16 FAR_FIELD
 17 WALL_TEMPERATURE
 18 TOTAL_TEMPERATURE
 19 TOTAL_PRESSURE
'''

GenerateLayout='''#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <array>

namespace patch
{
/* patch to make to_string work */
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}
}

using namespace std;
struct coord
{
    long double x;
    long double y;
    long double z;
    bool operator==(const coord& a) const
    {
          return (x == a.x && y == a.y && z == a.z);
    }
};

struct medians
{
    coord left;
    coord right;
    coord bottom;
    coord top;
    coord back;
    coord front;
};

struct vertices
{
    coord a;
    coord b;
    coord c;
    coord d;
    coord e;
    coord f;
    coord g;
    coord h;
};

struct faces
{
    std::string left,right,bottom,top,back,front;
};



class read_compute
{
public:
    int tot_grids;
    vector <int> imax, jmax, kmax;
    std::vector <std::string> grid_names;
    std::vector <vertices> all_vertices;
    std::vector <faces> all_faces;
    std::vector<std::array<int, 6>> mpi_class;
    void read_grid_names(std::ostream& block);//read file_names to get grids
    void compute_grid_medians(std::string str, std::ostream& block );
    void compute_interfaces();
    void write_to_file();
    void write_to_mapfile(std::ostream& map, int b1,int f1,int s11,int s12,int e11,
    int e12,int b2,int f2,int s21,int s22,int e21,int e22, bool dir_swap);
    void check_interface(int b1, int b2, int fb, int fc, 
        coord v11, coord v12, coord v13, coord v14,
        coord v21, coord v22, coord v23, coord v24, std::ostream& map);
    std::string layout_comment(std::string str);
};

void read_compute::write_to_file()
{

    string str = "layout.md",comnt_str;
    char *c = &str[0u];
    std::fstream output_file (c,fstream::out);
    if(!output_file.is_open())
    {
        cout << "Error: Couldn't open output file" << endl;
        return;
    }
    comnt_str = "BLOCK LAYOUT FILE";
    output_file << layout_comment(comnt_str);
    comnt_str = "==========================";
    output_file << layout_comment(comnt_str);
    comnt_str = "NUMBER OF PROCESSES";
    output_file << layout_comment(comnt_str);
    output_file << tot_grids << "\\n";
    comnt_str = "NUMBER OF ENTRIES PER PROCESS";
    output_file << layout_comment(comnt_str);
    output_file << "9\\n";
    comnt_str = "PROCESS_NO GRID BC_FILE IMIN IMAX JMIN JMAX KMIN KMAX";
    output_file << layout_comment(comnt_str);
    comnt_str = "===================================";
    output_file << layout_comment(comnt_str);
    for (int i =0; i< tot_grids; i++)
    {
    comnt_str = "PROCESS "+ patch::to_string(i);
    output_file << layout_comment(comnt_str);
    output_file <<setfill('0') <<setw(2) <<i << "  " << grid_names[i].substr(13,-1)<<"  "<<"bc_" << setw(2) << setfill('0') << i<<".md" <<"  ";
    output_file << setw(4) << all_faces[i].left << "  " << setw(4) << all_faces[i].right << "  "<< setw(4) << all_faces[i].bottom << "  ";
    output_file << setw(4) << all_faces[i].top << "  "<< setw(4) << all_faces[i].front << "  "<< setw(4) << all_faces[i].back << "  ";
    output_file << "\\n";
    }
    output_file.close();

}

void read_compute::compute_interfaces()
{

    bool dir_swap=false;
    vertices buf1,buf2;
    coord b_Va,b_Vb,b_Vc,b_Vd;
    coord b_Ve,b_Vf,b_Vg,b_Vh;
    coord c_Va,c_Vb,c_Vc,c_Vd;
    coord c_Ve,c_Vf,c_Vg,c_Vh;
    std::ofstream map;
    map.open("mapping.txt");
    map << "#   B1    F1    S1    E1    S2    E2    B2    F2    S1    E1    S2    E2    dir_swap" << endl;
    for(int i =0; i< tot_grids ; i++)
    {
        /*
        imin face a -> c -> e -> g
        imax face b -> d -> f -> h
        jmin face a -> b -> e -> f
        jmax face c -> d -> g -> h
        kmin face a -> b -> c -> d
        kmax face e -> f -> g -> h
        */
        // left face check
//        if(all_faces[i].left == "USER")
//        {
           //cout << "enterd 1";
            buf1 = all_vertices[i];
            for(int j = 0; j< tot_grids; j++)
            {
              if(i!=j){
                buf2 = all_vertices[j];
                b_Va = buf1.a;
                b_Vb = buf1.b;
                b_Vc = buf1.c;
                b_Vd = buf1.d;
                b_Ve = buf1.e;
                b_Vf = buf1.f;
                b_Vg = buf1.g;
                b_Vh = buf1.h;
                c_Va = buf2.a;
                c_Vb = buf2.b;
                c_Vc = buf2.c;
                c_Vd = buf2.d;
                c_Ve = buf2.e;
                c_Vf = buf2.f;
                c_Vg = buf2.g;
                c_Vh = buf2.h;
                //cout << c_buf1.x << " " << c_buf2.x << endl;
                check_interface(i,j,1,1,b_Va,b_Vc,b_Vg,b_Ve,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,1,2,b_Va,b_Vc,b_Vg,b_Ve,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,1,3,b_Va,b_Vc,b_Vg,b_Ve,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,1,4,b_Va,b_Vc,b_Vg,b_Ve,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,1,5,b_Va,b_Vc,b_Vg,b_Ve,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,1,6,b_Va,b_Vc,b_Vg,b_Ve,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,2,1,b_Vb,b_Vd,b_Vh,b_Vf,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,2,2,b_Vb,b_Vd,b_Vh,b_Vf,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,2,3,b_Vb,b_Vd,b_Vh,b_Vf,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,2,4,b_Vb,b_Vd,b_Vh,b_Vf,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,2,5,b_Vb,b_Vd,b_Vh,b_Vf,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,2,6,b_Vb,b_Vd,b_Vh,b_Vf,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,3,1,b_Va,b_Vb,b_Vf,b_Ve,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,3,2,b_Va,b_Vb,b_Vf,b_Ve,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,3,3,b_Va,b_Vb,b_Vf,b_Ve,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,3,4,b_Va,b_Vb,b_Vf,b_Ve,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,3,5,b_Va,b_Vb,b_Vf,b_Ve,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,3,6,b_Va,b_Vb,b_Vf,b_Ve,c_Ve,c_Vf,c_Vh,c_Vg,map);
                
                check_interface(i,j,4,1,b_Vc,b_Vd,b_Vh,b_Vg,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,4,2,b_Vc,b_Vd,b_Vh,b_Vg,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,4,3,b_Vc,b_Vd,b_Vh,b_Vg,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,4,4,b_Vc,b_Vd,b_Vh,b_Vg,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,4,5,b_Vc,b_Vd,b_Vh,b_Vg,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,4,6,b_Vc,b_Vd,b_Vh,b_Vg,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,5,1,b_Va,b_Vb,b_Vd,b_Vc,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,5,2,b_Va,b_Vb,b_Vd,b_Vc,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,5,3,b_Va,b_Vb,b_Vd,b_Vc,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,5,4,b_Va,b_Vb,b_Vd,b_Vc,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,5,5,b_Va,b_Vb,b_Vd,b_Vc,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,5,6,b_Va,b_Vb,b_Vd,b_Vc,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,6,1,b_Ve,b_Vf,b_Vh,b_Vg,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,6,2,b_Ve,b_Vf,b_Vh,b_Vg,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,6,3,b_Ve,b_Vf,b_Vh,b_Vg,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,6,4,b_Ve,b_Vf,b_Vh,b_Vg,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,6,5,b_Ve,b_Vf,b_Vh,b_Vg,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,6,6,b_Ve,b_Vf,b_Vh,b_Vg,c_Ve,c_Vf,c_Vh,c_Vg,map);

//                if(b_Va==c_Vb && b_Vc==c_Vd && b_Ve==c_Vf && b_Vg==c_Vh){
//                  all_faces[i].left =patch::to_string(j);
//                  all_faces[j].right =patch::to_string(i);
//                  swap=false;
//                  write_to_mapfile(map,i,1,1,jmax[i],1,kmax[i],j,3,1,jmax[j],1,kmax[j],swap);
//                  break;
//                }
//                if(b_Va==c_Vd && b_Vc==c_Vh && b_Ve==c_Va && b_Vg==c_Vf){
//                  all_faces[i].left =patch::to_string(j);
//                  all_faces[j].right =patch::to_string(i);
//                  swap=true;
//                  write_to_mapfile(map,i,1,1,jmax[i],1,kmax[i],j,3,1,jmax[j],kmax[j],1,swap);
//                  break;
//                }
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//
//                    continue;
//                }
//                c_buf1 = buf1.c;
//                c_buf2 = buf2.d;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                //cout << "enterd 10";
//                    continue;
//                }
//                c_buf1 = buf1.e;
//                c_buf2 = buf2.f;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.g;
//                c_buf2 = buf2.h;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                all_faces[i].left =patch::to_string(j);
//                all_faces[j].right =patch::to_string(i);
//                break;
              }
 //           }
        }

//        // top face check
//        if(all_faces[i].bottom == "USER")
//        {
//            buf1 = all_vertices[i];
//            for(int j = 0; j< tot_grids; j++)
//            {
//                buf2 = all_vertices[j];
//                c_buf1 = buf1.a;
//                c_buf2 = buf2.c;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.b;
//                c_buf2 = buf2.d;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.e;
//                c_buf2 = buf2.g;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.f;
//                c_buf2 = buf2.h;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                all_faces[i].bottom =patch::to_string(j);
//                all_faces[j].top =patch::to_string(i);
//                break;
//            }
//        }
//
//        // front face check
//        if(all_faces[i].front == "USER")
//        {
//            buf1 = all_vertices[i];
//            for(int j = 0; j< tot_grids; j++)
//            {
//                buf2 = all_vertices[j];
//                c_buf1 = buf1.e;
//                c_buf2 = buf2.a;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.g;
//                c_buf2 = buf2.c;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.f;
//                c_buf2 = buf2.b;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.h;
//                c_buf2 = buf2.d;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                all_faces[i].front =patch::to_string(j);
//                all_faces[j].back =patch::to_string(i);
//                break;
//            }
//        }
//
    }
    map.close();

}
void read_compute::read_grid_names(std::ostream& block)
{
    std::string str_buf;
    int count =0;
    std::ifstream file_buf;
    file_buf.open("grid_names");
    while(file_buf >> str_buf)
    {
        grid_names.push_back(str_buf);
        count++;
    }
    tot_grids = count;
    block << "MBLK \t"<<tot_grids << endl;
}

void read_compute::compute_grid_medians(std::string str, std::ostream& block)
{
    char *cr = &str[0u];
    std::ifstream grid;
    grid.open(cr);
    cout << cr << endl;
    int imx =51,jmx =51, kmx =51;
    coord a,b,c,d,e,f,g,h,buf;
    medians buf_med;
    grid >> imx >> jmx >> kmx;
    imax.push_back(imx);
    jmax.push_back(jmx);
    kmax.push_back(kmx);
    block << imx << "\t" << jmx << "\t" << kmx << endl;
    grid >> a.x >> a.y >> a.z;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> b.x >> b.y >> b.z;
    for(int j = 2; j < jmx; j++)
    {
        for(int i =1 ; i<=imx ; i++)
        {
            grid >> buf.x >> buf.y >> buf.z;
        }
    }
    grid >> c.x >> c.y >> c.z;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> d.x >> d.y >> d.z;

    //skip in k direction
    for (int k =2; k< kmx; k++)
        for(int j = 1; j <= jmx; j++)
        {
            for(int i =1 ; i<=imx ; i++)
            {
                grid >> buf.x >> buf.y >> buf.z;
            }
        }

    grid >> e.x >> e.y >> e.z;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> f.x >> f.y >> f.z;
    for(int j = 2; j < jmx; j++)
    {
        for(int i =1 ; i<=imx ; i++)
        {
            grid >> buf.x >> buf.y >> buf.z;
        }
    }
    grid >> g.x >> g.y >> g.z;
    //cout << g.x << " " << g.y << " "<< g.z << endl;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> h.x >> h.y >> h.z;

    // store vertices
    vertices buf_vert;
    buf_vert.a = a;
    buf_vert.b = b;
    buf_vert.c = c;
    buf_vert.d = d;
    buf_vert.e = e;
    buf_vert.f = f;
    buf_vert.g = g;
    buf_vert.h = h;
    all_vertices.push_back(buf_vert);
    // store faces
    faces buf_face;
    buf_face.back = "USER";
    buf_face.front = "USER";
    buf_face.left = "USER";
    buf_face.right = "USER";
    buf_face.top = "USER";
    buf_face.bottom = "USER";
    all_faces.push_back(buf_face);
    //j10 : faces types for mpi sequence
    std::array<int, 6> buf_mpi_class={-1,-1,-1,-1,-1,-1};
    mpi_class.push_back(buf_mpi_class);
    grid.close();
}
std::string read_compute::layout_comment(string str)
{
return "## " + str + "\\n";
}

void read_compute::write_to_mapfile(std::ostream& map, int b1,int f1,int s11,int s12,int e11,
    int e12,int b2,int f2,int s21,int s22,int e21,int e22, bool dir_swap){
map<<"  " << setw(4) << b1 <<"  "<<setw(4) << f1  <<"  "<<setw(4) << s11<<"  "<<setw(4) << s12 
   <<"  " << setw(4) << e11<<"  "<<setw(4) << e12 <<"  "<<setw(4) << b2 <<"  "<<setw(4) << f2
   <<"  " << setw(4) << s21<<"  "<<setw(4) << s22 <<"  "<<setw(4) << e21<<"  "<<setw(4) << e22
   <<"    "<<dir_swap << "  " << mpi_class[b1][f1]<<endl;
}

void read_compute::check_interface(int i, int j, int f1, int f2, coord v11, coord v12, coord v13, coord v14,
        coord v21, coord v22, coord v23, coord v24, std::ostream& map){
  bool dir_swap=false;
  bool found=false;
  int  type=1;
  int  s11=1,s12=1,s21=1,s22=1;
  int  e11,e12,e21,e22;
  if(v11==v21 && v12==v22 && v13==v23 && v14==v24){
    dir_swap=false; found=true; type=1;
  }
  if(v11==v22 && v12==v23 && v13==v24 && v14==v21){
    dir_swap=true; found=true; type=2;
  }
  if(v11==v22 && v12==v21 && v13==v24 && v14==v23){
    dir_swap=true; found=true; type=2;
  }
  if(v11==v23 && v12==v24 && v13==v21 && v14==v22){
    dir_swap=false; found=true; type=3;
  }
  if(v11==v24 && v12==v21 && v13==v22 && v14==v23){
    dir_swap=true; found=true; type=4;
  }
//  if(v11==v24 && v12==v23 && v13==v22 && v14==v21){
//    dir_swap=false; found=true; type=4;
  if(found==true){
    if(f1==1){ all_faces[i].left   =patch::to_string(j); e11=jmax[i]; e12=kmax[i];}
    if(f1==2){ all_faces[i].right  =patch::to_string(j); e11=jmax[i]; e12=kmax[i];}
    if(f1==3){ all_faces[i].bottom =patch::to_string(j); e11=imax[i]; e12=kmax[i];}
    if(f1==4){ all_faces[i].top    =patch::to_string(j); e11=imax[i]; e12=kmax[i];}
    if(f1==5){ all_faces[i].front  =patch::to_string(j); e11=imax[i]; e12=jmax[i];}
    if(f1==6){ all_faces[i].back   =patch::to_string(j); e11=imax[i]; e12=jmax[i];}
    if(f2==1){ all_faces[j].left   =patch::to_string(i); e21=jmax[j]; e22=kmax[j];}
    if(f2==2){ all_faces[j].right  =patch::to_string(i); e21=jmax[j]; e22=kmax[j];}
    if(f2==3){ all_faces[j].bottom =patch::to_string(i); e21=imax[j]; e22=kmax[j];}
    if(f2==4){ all_faces[j].top    =patch::to_string(i); e21=imax[j]; e22=kmax[j];}
    if(f2==5){ all_faces[j].front  =patch::to_string(i); e21=imax[j]; e22=jmax[j];}
    if(f2==6){ all_faces[j].back   =patch::to_string(i); e21=imax[j]; e22=jmax[j];}

    if(type==2){swap(s21,e21);type=1;}
    if(type==3){swap(s21,e21); swap(s22,e22);}
    if(type==4){swap(s22,e22);}
    if(mpi_class[i][f1]<0 && mpi_class[j][f2]<0){
      mpi_class[i][f1]=0;
      mpi_class[j][f2]=1;
    }
    write_to_mapfile(map,i,f1,s11,e11,s12,e12,j,f2,s21,e21,s22,e22,dir_swap);
  }
  
}

int main()
{
    std::ofstream block;
    block.open("blocking.dat");
    read_compute handler;
    handler.read_grid_names(block);
    block << "imx \t jmx \t kmx \\n";
    for(int i =0; i< handler.tot_grids; i++)
    {
        cout << "Wait!! Reading grid "<< i+1 << endl;
        handler.compute_grid_medians(handler.grid_names[i], block);
    }
    cout << "computing adjacent blocks" << endl;
    handler.compute_interfaces();
    handler.write_to_file();
    block.close();
    //cout << handler.all_faces[0].left << handler.all_faces[0].right << endl;
    //cout << handler.all_vertices.size();
    //cout << handler.grid_names[1] << handler.tot_grids << endl;
    return 0;
}
'''

GenerateBC='''#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <utility>

namespace patch
{
/* patch to make to_string work */
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}
}

using namespace std;

struct all_bc
{
    int imin;
    int imax;
    int jmin;
    int jmax;
    int kmin;
    int kmax;

};

class read_write
{
public:
    std::map <int,std::string> bc_values;
    std::map <int,std::string> bc_list;
    void read_bc_values();
    void read_bc_list();
    void read_layout();
    void write_bc_file(std::string, all_bc);
    int check_for_periodic_bc(std::ostream& periodic, int bc,  int count, int face );
    std::string parse_bc(std::string);

};
std::string read_write::parse_bc(std::string bc)
{
    std::string buf ="",net_string="";
    int bc_num;
    for(int i =0; i< bc.length(); i++)
    {
        if(bc[i] == '-')
        {

            bc_num =atoi(buf.c_str());
            net_string =net_string+ "- " + bc_values[bc_num] +"\\n";
            buf = "";
        }
        else
            buf = buf+bc[i];
    }
    bc_num =atoi(buf.c_str());
    net_string =net_string+ "- " + bc_values[bc_num] +"\\n\\n";
    buf = "";
    return net_string;
}

void read_write::write_bc_file(std::string bc_file, all_bc buf)
{
    string str = "../bc/"+bc_file;
    char *c = &str[0u];
    std::fstream output_file (c,fstream::out);
    if(!output_file.is_open())
    {
        cout << "Error: Couldn't open output file" << endl;
        return;
    }
    output_file << "BOUNDARY CONDITIONS CONFIGURATION\\n=================================\\n\\n";

    std::map<int,std::string>::iterator it;
    std::string bc;
    int bc_num;

    // imin
    bc_num = buf.imin;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# imn\\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

     // imax
    bc_num = buf.imax;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# imx\\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

         // jmin
    bc_num = buf.jmin;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# jmn\\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

         // jmax
    bc_num = buf.jmax;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# jmx\\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

             // kmin
    bc_num = buf.kmin;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# kmn\\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

         // kmax
    bc_num = buf.kmax;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# kmx\\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);
    output_file << "FIN\\n";



}

void read_write::read_layout()
{
    std::ofstream periodic;
    periodic.open("periodic.txt");
    periodic << "   B1   B2   F1   F2 "<<endl;
    fstream fin("layout.md");
    std::string line;
    int num,count=0,tot_proc,tot_ent;
    int proc_id;
    std::string imin,imax,jmin,jmax,kmin,kmax;
    //std::string imin="USER",imax="USER",jmin="USER",jmax="USER",kmin="USER",kmax="USER";
    std::string grid,bc_file;
    all_bc buf;
    while(getline(fin, line))
    {
        //the following line trims white space from the beginning of the string
        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));

        if(line[0] == '#') continue;//ignore lines starting with #
        count++;

        if(count > 2)
        {
            // boundary condition line
            stringstream(line) >> proc_id >> grid >> bc_file >> imin >> imax >> jmin >> jmax >> kmin >> kmax;
             if(imin == "USER" || imax == "USER" || jmin =="USER" ||jmax =="USER" || kmin == "USER" || kmax =="USER")
            {
               cout << "ERROR: ONE OR MORE USER BCs ARE YET TO BE GIVEN FOR PROCESS " << proc_id << endl;
               cout << "EROOR: NOT GENERATING BC FILE " << bc_file << endl;
               continue;
            }

            buf.imin = atoi(imin.c_str());
            buf.imax = atoi(imax.c_str());
            buf.jmin = atoi(jmin.c_str());
            buf.jmax = atoi(jmax.c_str());
            buf.kmin = atoi(kmin.c_str());
            buf.kmax = atoi(kmax.c_str());
            buf.imin = check_for_periodic_bc(periodic, buf.imin, count-3, 1);
            buf.jmin = check_for_periodic_bc(periodic, buf.jmin, count-3, 3);
            buf.kmin = check_for_periodic_bc(periodic, buf.kmin, count-3, 5);
            buf.imax = check_for_periodic_bc(periodic, buf.imax, count-3, 2);
            buf.jmax = check_for_periodic_bc(periodic, buf.jmax, count-3, 4);
            buf.kmax = check_for_periodic_bc(periodic, buf.kmax, count-3, 6);
            write_bc_file(bc_file,buf);

        }
        else if(count == 1)
            stringstream(line) >> tot_proc;
        else if(count == 2)
            stringstream(line) >> tot_ent;
    }
    //cout << tot_proc << " " << tot_ent <<endl;

}

int read_write::check_for_periodic_bc(std::ostream& periodic, int bc, int count, int face){
  // get the first digit
  std::string s = std::to_string(bc);
  int periodicF2;
  if(bc<0 && s.length()>2 && s.substr(1,1)==std::to_string(2)) {
    //write periodic bc control file
    if(face == 1) {periodicF2 = 2;}
    if(face == 2) {periodicF2 = 1;}
    if(face == 3) {periodicF2 = 4;}
    if(face == 4) {periodicF2 = 3;}
    if(face == 5) {periodicF2 = 6;}
    if(face == 6) {periodicF2 = 5;}
    periodic << " " << setw(4) << count << " " << setw(4) << s.substr(2) 
             << " " << setw(4) << face << " " << setw(4) << periodicF2 << endl;
    return -10;
    }
  else{ 
    return bc;
    }
  }

void read_write::read_bc_values()
{
    fstream fin("bc_values");
    std::string line,bc;
    int num;
    while(getline(fin, line))
    {
        //the following line trims white space from the beginning of the string
        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));

        if(line[0] == '#') continue;//ignore lines starting with #


        stringstream(line) >> num >> bc;
        bc_values[num] = bc;
        //cout << "Data: " << num << bc << endl;
    }
}

void read_write::read_bc_list()
{
    fstream fin("bc_list");
    std::string line,bc;
    int num;
    while(getline(fin, line))
    {
        //the following line trims white space from the beginning of the string
        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));

        if(line[0] == '#') continue;//ignore lines starting with #


        stringstream(line) >> num >> bc;
        bc_list[num] = bc;
        //cout << "Data: " << num << bc << endl;
    }
}


int main()
{
    read_write handler;
    handler.read_bc_values();
    handler.read_bc_list();
    handler.read_layout();
    return 0;
}
'''

BashRun='''#!/bin/bash

ex='bin/FEST3D'
runlog='time_directories/aux/out'
total_process='''+str(NumberOfBlocks)+''' #"$2" #parallel running - number of process

echo "Run: `date +%Y/%m/%d-%H/%M/%S`" | tee -a $runlog

if [ -f $ex ]; then
    echo 'Running solver.' | tee -a $runlog
    echo >> $runlog
    echo 'Solver output:' >> $runlog
    echo >> $runlog
    mpiexec.hydra -np $total_process ./$ex >> $runlog
fi

echo | tee -a $runlog
echo 'End of run script.' | tee -a $runlog
 '''


Compiler='''g++ -std=c++11 $1
'''
NameFiller='''#!/bin/bash
total_process="$1"
echo -n > grid_names
for((i =0; i<total_process;i++ ))
do
printf -v j "%02d" $i
echo "../gridfiles/grid_$j.txt" >> grid_names
done
'''
GnuplotPlot='''set term x11 1 noraise
set logscale y
set grid
#set yrange [1e-8:1e3]
plot for [col=2:'''+str(len(ResidualControl['Out'])+1)+'''] 'time_directories/aux/resnorm' using 1:col with lines title columnheader
pause 5; unset output; refresh; reread;
'''
VTKInputFiller='''#!/bin/bash
total_process="$1"
Time_stamp="$2"
printf -v j "%04d" $2
file="time_directories/$j/list.visit"
echo -n > $file
echo !NBLOCKS $total_process >> $file
for((i =0; i<total_process;i++ ))
do
printf -v j "%02d" $i
echo "process_$j.vtk" >> $file
done
'''


def SetBC(RootDir,BoundaryConditions):
    OldPath = RootDir+"/system/mesh/layout/layout.md"
    NewPath = RootDir+"/system/mesh/layout/layout.tp"
    with open(NewPath, "w+") as NewFile:
        with open(OldPath, "r") as OldFile:
            for line in OldFile:
                if(len(line)>5 and not(line[0]=='#')):
                    data = line.split()
                    domain = data[3:]
                    for num, boundary in enumerate(domain):
                        if(boundary=='USER'):
                            domain[num] = "{0:-04d}".format(int(BoundaryConditions[num]))
                    data[3:] = domain
                    CorrectedLine = "  ".join(data)+'\n'
                else:
                    CorrectedLine = line
                NewFile.write(CorrectedLine)
    os.remove(OldPath)
    shutil.move(NewPath, OldPath)

WriteFullFile(RunDir, BCListFilename, BCList)
WriteFullFile(RunDir, BCValueFilename, BCValue)
WriteFullFile(RunDir, GenLayoutFilename, GenerateLayout)
WriteFullFile(RunDir, GenBCFilename, GenerateBC)
WriteFullFile(RunDir, CompilerFilename, Compiler)
WriteFullFile(RunDir, NameFillerFilename, NameFiller)
WriteFullFile(RunDir, BashRunFilename, BashRun)
WriteFullFile(RunDir, GnuplotFilename, GnuplotPlot)
if Control['OutputFileFormat'] == 'vtk':
  WriteFullFile(RunDir, VTKFillerFilename, VTKInputFiller)

# copy Gridfile to particular folder
GridfilesFolder=RunDir+"/system/mesh/gridfiles/"
files = os.listdir(GridDir)
for f in files:
    shutil.copy2(GridDir+'/'+f, GridfilesFolder)

LayoutDir=RunDir+"/system/mesh/layout/"
p = subprocess.Popen(["bash", FillerFile, str(NumberOfBlocks)], cwd=LayoutDir)
p = subprocess.Popen(["bash", CompileFile, GenLayoutFile], cwd=LayoutDir)
p.wait()
p = subprocess.Popen(["./a.out"], cwd=LayoutDir)
p.wait()
SetBC(RunDir, BoundaryConditions)
p = subprocess.Popen(["bash", CompileFile, GenBCFile], cwd=LayoutDir)
p.wait()
p = subprocess.Popen(["./a.out"], cwd=LayoutDir)
p.wait()
p = subprocess.Popen(["ln", "-fs", os.environ['FEST3D'], RunDir+"/bin/FEST3D"])
p.wait()
