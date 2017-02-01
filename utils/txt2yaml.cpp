/* To be improved.
* Expected syntax to generate the correct yaml file from command line:
* ./generator 'inputFile' 'outputFile'
* COMPILE WITH: c++ -o yaml_datafiles_generator.cpp yaml_datafiles_generator.exe
*********************************************************************************/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>

using namespace std;

int main(int argc, char* argv[])
{
    
    if(argc%2==0){
        cout << "Error: wrong number of arguments!\n";
        return 1;
    }
    
    else{
        cout << "argc  = "<<argc << "argv[0]" << argv[0]<<endl;
        int a;
        for(a=0; a<(int) (argc-1)/2; a++){
    
            string ind_var, dep_var, err, sys;
            string ind_var_unit, dep_var_unit, err_unit;
            string line;

            ifstream dummy(argv[2*a+1]);

            int n=0;
cout<<"reading"<<endl;
            while(!dummy.eof()){
                getline(dummy, line);
                n=n+1;
            }
            n=n-3;

            dummy.close();
cout<<"done, n = "<<n << endl;
cout << "read "<< argv[2*a+1]<< " translate to "<<argv[2*a+2]<< endl;
            ifstream data(argv[2*a+1]);
            ofstream yaml(argv[2*a+2]);

            data >> ind_var >> dep_var >> err >> sys;
            data >> ind_var_unit >> dep_var_unit >> err_unit >> err_unit;

            double x[n], y[n], e[n], s[n];

            int i=0;
cout<<"2 reading"<<endl;
            while(!data.eof()){
            //for(int i = 2; i < n; i++){
                data >> x[i] >> y[i] >> e[i] >> s[i];
                //cout << i << x[i] << y[i] << e[i];
                i=i+1;
            }
cout<<"done"<<endl;
            yaml << "name: '" << argv[2*a+2] << "'" << "\n";
            yaml << "independent_variables:" << "\n";
            yaml << " - header: {name: " << ind_var << ", units: " << ind_var_unit << "}" << "\n";
            yaml << "   values:" << "\n";

            for(i=0; i<n; i++){
                yaml << "    - value : " << x[i] << "\n";
            }
         cout<< "here"<<endl;       
            yaml << "dependent_variables:" << "\n";
            yaml << " - header: {name: " << dep_var << ", units: " << dep_var_unit << "}" << "\n";
            yaml << "   qualifiers:" << "\n";
            yaml << "    - {name: '', value: ''}" << "\n";
            yaml << "   values:" << "\n";
                
            for(i=0; i<n; i++){
                yaml << "     - value: " << y[i] << "\n";
                yaml << "       errors:" << "\n";
                yaml << "         - {symerror: " << e[i] << ", label: stat}" << "\n";
                yaml << "         - {symerror: " << s[i] << ", label: sys}" << "\n";
            }

            data.close();
            yaml.close();
        }
    }
    
	return 0;
}
