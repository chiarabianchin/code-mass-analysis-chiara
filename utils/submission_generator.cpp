/* To be improved.
* Expected syntax to generate the correct yaml file from command line:
* ./generator 'inputFile' 'outputFile'
* COMPILE WITH: c++ -o yaml_submissionfile_generator.cpp yaml_submissionfile_generator.exe
*********************************************************************************/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>

using namespace std;

int main(int argc, char* argv[])
{
	//ppb
    string ppbreaction  = "P PB --> JET X";
    double ppbcme = 5020;
	string pbpbreaction = "PB PB --> JET X";
	string centppb = "MB";
	double pbpbcme = 2760;
	//pbpb
	double radius = 0.4;
	string centpbpb = "0.0-10.0";
	
	string obsdist = "DN/DLIGHT-JET-MASS", obsmean = "MEANLIGHT-JET-MASS", obs = obsdist;
		
	//observables
	string var = "DN/DLIGHT-JET-MASS", ratio = "RATIO DN/DLIGHT-JET-MASS";
	
	ifstream table_files(argv[1]);
    ofstream submission("submission.yaml");
    
    submission << "comment: | #Jet mass in p--Pb and Pb--Pb collisions measured by ALICE" << "\n";
    submission << "\n\n\n\n\n";
    
    int i=1;
    string filename;
    
    string reaction, cent;
    double cme;
    
    while(1){
    	bool isppb = false, ispbpb = false, isratio = false;
    	
        table_files >> filename;
        cout << filename << endl;
        if(table_files.eof()) break;
        
        submission << "---\n\n";
        submission << "#This is Table " << i << "\n";
        submission << "name: 'Table " << i << "'" << "\n";
        submission << "location: data from Fig. " << "\n";
        submission << "description: description to be inserted" << "\n";
        submission << "keywords: " << "\n";
        
        if(filename.find("pPb")  != std::string::npos) {
        	cout << "pPb collisions"<<endl;
        	isppb = true;
        	reaction = ppbreaction;
        	cent = centppb;
        	cme = ppbcme;
        }
        if(filename.find("PbPb") != std::string::npos) {
        	cout << "PbPb collisions"<<endl;
        	ispbpb = true;
        	reaction = pbpbreaction;
        	cent = centpbpb;
        	cme = pbpbcme;
        }
        if(isppb && ispbpb) isratio = true;
        if(filename.find("Mean") != std::string::npos){
        	obs = obsmean;
        }
        if(isratio) {
        submission << " - { name: reactions, values: [" << ppbreaction << " and "<< pbpbreaction << "]}" << "\n";
        submission << " - { name: observables, values: [" << ratio << "]}" << "\n";
        submission << " - { name: cmenergies, values: [" << ppbcme << " and "<< pbpbcme << "]}" << "\n";
        } else {
        	submission << " - { name: reactions, values: [" << reaction << "]}" << "\n";
        	submission << " - { name: observables, values: [" << obs << "]}" << "\n";
        	submission << " - { name: cmenergies, values: [" << cme << "]}" << "\n";
        }
        submission << "- {name: phrases, values: ['Jet mass']"<< endl;
        submission << "data_file: " << filename << "\n\n\n";
        
        i=i+1;
    }

    return 0;
}
