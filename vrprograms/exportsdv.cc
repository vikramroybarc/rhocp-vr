// This Program is Written to Create a Input FIle for exporting the SDV Data From Moose Simulations
#include<iostream>
#include<cstring>
#include<fstream>

void Usage(std::string);

int main(int argc, char* argv[]){
    int n,  arg_no;
    int pflag{0}, nflag{0}, oflag{0}, fflag{0};
    std::string filename;
    

if (argc == 1){
    Usage(argv[0]);
    exit(0);
}

// Get Command Line Params
for (int i = 1; i < argc; i++){
    if (!strcmp(argv[i], "-h")){
        Usage(argv[0]);
        exit(0);
    }
    else if (!strcmp(argv[i], "-n")){
       n = std::stoi(argv[++i]);
       arg_no = i+1;
       nflag = 1; 
       break;
    }
    else {
        Usage(argv[0]);
        exit(1);
    }      
}

std::string prefix[n];
int sdv_start[n], sdv_end[n], sdvflag{0};
std::string order[n];



for (int i = arg_no; i < argc; i++){
    if (!strcmp(argv[i], "-pre")){
        for(int j = 0; j < n; j++){
            prefix[j] = argv[++i];            
        }
        pflag = 1;
    }
    else if (!strcmp(argv[i], "-s")){
        for(int j = 0; j < n; j++){
            sdv_start[j] = std::stoi(argv[++i]);
            sdv_end[j] = std::stoi(argv[++i]);
        }
        sdvflag = 1;
    }
    else if (!strcmp(argv[i], "-ord")){
        for(int j = 0; j < n; j++){
            order[j] = argv[++i];
        }
        oflag = 1;
    }  
    else if (!strcmp(argv[i], "-f")){
        filename = argv[++i];
        fflag = 1;
    }
    else {
        Usage(argv[0]);
        exit(1);
    }          
    }

    if (!fflag){
        filename = "moose_sdv.i";
        fflag = 1;
    }

    // Open OutPut File for Writing
    std::fstream ofile(filename, std::ios::out);
    ofile << "[AuxVariables]" << std::endl;
    int ntype[n]{0};

    for (int i = 0; i < n; i++){
        int typen_no = sdv_end[i] - sdv_start[i] +1 ;
        
        std::string varorder;
        if((order[i]) == "C"){
            varorder = "CONSTANT";
        }
        else if (order[i] == "F"){
            varorder = "FIRST";
        }
        else {
            varorder = "FIRST";
        }

        std::string varfamily{"MONOMIAL"};

        // Define AUXVars
        for(int j = 0; j < typen_no; j++){
            std::string varname{prefix[i]};
            varname.append("_" + std::to_string(j+1));
            ofile << "\t [" + varname + "]\n";
            ofile << "\t \t order = " + varorder + "\n"; 
            ofile << "\t \t family = " + varfamily + "\n"; 
            ofile << "\t []\n" << std::endl;

        }
         
    }
    ofile << "[]" << std::endl;
    
    // Define AUXKernals
    ofile << "\n\n[AuxKernels]" << std::endl;
    for (int i = 0; i < n ; i++){
    int typen_no = sdv_end[i] - sdv_start[i] + 1;

        
        for(int j = 0; j < typen_no; j++){
            std::string varname{prefix[i]};
            varname.append("_" + std::to_string(j+1));
            ofile << "\t [" + varname + "]\n";
            ofile << "\t \t  type = StateVariable \n"; 
            ofile << "\t \t variable = " + varname + "\n";
            ofile << "\t \t sdv_id = " + std::to_string(sdv_start[i]+j) << std:: endl;
            ofile << "\t \t execute_on = timestep_end\n";
            ofile << "\t []" << std::endl;
        }

    }
    ofile << "[]" << std::endl;

    // Define PostProcessors
    ofile << "\n\n [Postprocessors]" << std::endl;
    for (int i = 0; i < n; i++){
        int typen_no = sdv_end[i] - sdv_start[i] + 1;

        for (int j = 0; j < typen_no ; j++){
            std::string varname{prefix[i]};
            varname.append("_" + std::to_string(j+1));
            ofile << "\t [" +  varname + "]\n";
            ofile << "\t \t type = ElementAverageValue\n";
            ofile << "\t \t variable = " + varname << std::endl;
            ofile << "\t []" << std::endl;
        }
    }
    ofile << "[]" << std::endl;

    return 0;


}

void Usage(std::string programname){
    std::cout << programname << std::endl;
    std::cout << "\t [-n <sdv types>]" << std::endl;
    std::cout << "\t [-pre <prefix_sdv_type_1> <prefix_sdv_type_2> ... <prefix_sdv_type_n>]" << std::endl;
    std::cout << "\t [-s <sdv_start_no_type_1> <sdv_end_no_type_1> ... <sdv_start_no_type_n> <sdv_end_no_type_n>]" << std::endl;
    std::cout << "\t [-ord <auxvar_order_type1> <auxvar_order_type2> <auxvar_order_typen>]" << std::endl;
    std::cout << "\t [-f <filename>]" << std::endl;
    std::cout << "\t [-h]" << std::endl;
    std::cout << "\t Where :" << std::endl;
    std::cout << "\t sdv_types = No of Different Types of SDVs like: " 
              << "\n\t ClimbStrain, GlideStrain etc" << std::endl;              
    std::cout << "\t <prefix_sdv_type_1> = Prefix in AuxVar Name for SDV type 1"
              << "\n\t for eg: for glidestrain prefix=glidestrain" << std::endl;
    std::cout << "\t <sdv_start_no_type_1> <sdv_end_no_type_1> available in source file of that Material" << std::endl;
    std::cout << "\t -ord = C for AuxVars of Const Order and Monomial Family" << std::endl;
    std::cout << "\t -ord = F for AuxVars of First Order and Monomial Family" << std::endl;
    std::cout << "\t filename = Name of File in which all data regarding auxvars and auxkernals is to written" << std::endl;
    std::cout << "\t -h Prints Out the Help Information" << std::endl;
    std::cout << "\t example syntax: \n";
    std::cout << "\t exportsdv -n 2 -pre gstrain cstrain -s 100 124 148 172 -ord F F" << std::endl;

}

