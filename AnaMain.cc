#include "DataAnalysis.h"

//compilation
//g++ `root-config --glibs --cflags ` DataAnalysis.o main.C 

int main(){

 DataAnalysis ana;
 cout<<"Data loaded"<<endl;
 ana.LoadData("angle.txt",2);
 cout<<"Data loaded"<<endl;
 ana.Fit("cos2","angle","Nb coinc",2);
 ana.Fit("cosn","angle","Nb coinc",2);
 ana.Fit("cos2Gaus","angle","Nb coinc",2);
 cout<<"Data loaded"<<endl;

 return 0;
}
