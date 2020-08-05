#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdio.h> 
using namespace std;
struct node{
    bool isLeaf;
    string attr;
    double splitval;
    node* left;
    node* right;
    string label;
    node(bool isLeaf,string attr,double splitval,string label){
        isLeaf = false;
        attr = "";
        splitval = 0;
        label = "";
    }
};
struct pro{
    double l = 0.0;;
    double r = 0.0;
};
struct dt{
        string bestAttr = "";
        double bestsplitval = 0;
};
typedef struct pro Struct;
Struct cal_pro(vector<double> small,vector<double> big,string credit,int n,vector<string> smallC,vector<string> bigC){
    Struct s;
    s.l= 0;
    s.r = 0;
    for(int i = 0; i < small.size();i ++){
        if(smallC.at(i) == credit){
            s.l ++;
        }
    }
    for(int i = 0; i < big.size();i ++){
        if(bigC.at(i) == credit){
            s.r ++;
        }
    }
    s.l = s.l/small.size();
    s.r = s.r/big.size();
    return s;
}
double isNan(double a){
    if(isnan(a)){
        a = 0;
        return a;
    }else{
        a = a;
        return a;
    }
}
int attr_pos(string n){
        if(n == "WC_TA"){
            return 0;
        }else if(n == "RE_TA"){
            return 1;
        }else if(n == "EBIT_TA"){
            return 2;
        }else if(n == "MVE_BVTD"){
            return 3;
        }else if(n == "S_TA"){
            return 4;
        }
}
double cal_E(string data[],int n){
    n = n -1;
    double AAA = 0;
    double AA = 0;
    double A = 0;
    double BBB = 0;
    double BB = 0;
    double B = 0;
    double C = 0;
    for(int i = 0; i < n; i++){
        if(data[i] == "AAA"){
            AAA ++;
        }else if(data[i] == "AA"){
            AA ++;
        }else if(data[i] == "A"){
            A ++;
        }else if(data[i] == "BBB"){
            BBB ++;
        }else if(data[i] == "BB"){
            BB ++;
        }else if(data[i] == "B"){
            B ++;
        }else if(data[i] == "CCC"){
            C ++;
        }
    }
    double result = -1*(isNan((AAA/n)*log2(AAA/n))+isNan((AA/n)*log2(AA/n))+isNan((A/n)*log2(A/n))+isNan((BBB/n)*log2(BBB/n))+isNan((BB/n)*log2(BB/n))+isNan((B/n)*log2(B/n))+isNan((C/n)*log2(C/n)));
    return result;
}
bool is_leaf(string data[][6],int n){
    int i = 1;
    if(n == 2){
        return true;
    }
    while(i < n-1){
        if(data[i][5] != data[i+1][5]){
            return false;
        }
        i++;
    }
    return true;
}
bool is_Null(string data[][6]){
    int i = 1;
    if(data[1][5] == ""){
        return true;
    }
    return false;
}
string cal_majority(string data[][6],int n){
    n = n -1;
    double AAA = 0;
    double AA = 0;
    double A = 0;
    double BBB = 0;
    double BB = 0;
    double B = 0;
    double C = 0;
    int dataI[7] = {0,0,0,0,0,0,0};
    string dataS[7] = {"AAA","AA","A","BBB","BB","B","CCC"};
    for(int i = 0; i < n; i++){
        if(data[i][5] == "AAA"){
            dataI[0] ++;
        }else if(data[i][5] == "AA"){
            dataI[1] ++;
        }else if(data[i][5] == "A"){
            dataI[2] ++;
        }else if(data[i][5] == "BBB"){
            dataI[3] ++;
        }else if(data[i][5] == "BB"){
            dataI[4] ++;
        }else if(data[i][5] == "B"){
            dataI[5] ++;
        }else if(data[i][5] == "CCC"){
            dataI[6] ++;
        }
    }
    string gg;
    gg = dataS[0];
    for(int i = 0; i < 6 ; i ++){
        if(dataI[0] < dataI[i+1]){
            gg = dataS[i+1];
            dataI[0] = dataI[i+1];
        }else{
            
        }
    }
    return gg;
}
dt choose_split(string dataTemp[][6],int n){
    double gain = 0;
    double bestgain = 0;
    dt ret;
    ret.bestAttr = "";
    double data[n-1][6];
    string data_credit[n-1];
    double splitval = 0.0;
    for(int i = 0; i < 5; i ++){
        for(int j = 0;j<n-1; j ++){
            data[j][i] = stod(dataTemp[j+1][i]);
            data_credit[j] = dataTemp[j+1][5];
        }
    }
    double E = cal_E(data_credit,n);
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < n-1; j ++){
            data_credit[j] = dataTemp[j+1][5];
        }
        for(int m = 0; m < n-2; m ++){
            for(int j = 0; j < n-2; j ++){
                if(data[j][i] > data[j+1][i]){
                    double temp = data[j+1][i];
                    data[j+1][i] = data[j][i];
                    data[j][i] = temp;
                    string t = data_credit[j+1];
                    data_credit[j+1] = data_credit[j];
                    data_credit[j] = t;                    
                }
            }
        }
        for(int j = 0; j < n-1; j ++){
            double gain = 0;
            splitval = 0.5*(data[j][i]+data[j+1][i]);
            int tmpLeft =0;
            int tmpRight = 0;
            vector<double> small;
            vector<double> big;
            vector<string> smallC;
            vector<string> bigC;
            int count = 0;
            for(int m = 0; m < n-1; m ++){
                if(data[m][i] <= splitval){
                    small.push_back(data[m][i]);
                    smallC.push_back(data_credit[m]);
                }else{
                    big.push_back(data[m][i]);
                    bigC.push_back(data_credit[m]);
                }
            }
            Struct AAA;
            AAA = cal_pro(small,big,"AAA",n,smallC,bigC);
            Struct AA;
            AA = cal_pro(small,big,"AA",n,smallC,bigC);
            Struct A;
            A = cal_pro(small,big,"A",n,smallC,bigC);
            Struct BBB;
            BBB = cal_pro(small,big,"BBB",n,smallC,bigC);
            Struct BB;
            BB = cal_pro(small,big,"BB",n,smallC,bigC);
            Struct B;
            B = cal_pro(small,big,"B",n,smallC,bigC);
            Struct C;
            C = cal_pro(small,big,"CCC",n,smallC,bigC);
            double I_left = -1*(isNan(AAA.l*log2(AAA.l)) + isNan(AA.l*log2(AA.l)) + isNan(A.l*log2(A.l)) + isNan(BBB.l*log2(BBB.l)) + isNan(BB.l*log2(BB.l)) + isNan(B.l*log2(B.l)) + isNan(C.l*log2(C.l)) );
        
            double I_right = -1*(isNan(AAA.r*log2(AAA.r)) + isNan(AA.r*log2(AA.r)) + isNan(A.r*log2(A.r)) + isNan(BBB.r*log2(BBB.r)) + isNan(BB.r*log2(BB.r)) + isNan(B.r*log2(B.r)) + isNan(C.r*log2(C.r)) );
            gain = E - ((AAA.l*small.size())/(n-1))*I_left - ((AA.l*small.size())/(n-1))*I_left - ((A.l*small.size())/(n-1))*I_left - ((BBB.l*small.size())/(n-1))*I_left - ((BB.l*small.size())/(n-1))*I_left - ((B.l*small.size())/(n-1))*I_left - ((C.l*small.size())/(n-1))*I_left - ((AAA.r*big.size())/(n-1))*I_right - ((AA.r*big.size())/(n-1))*I_right - ((A.r*big.size())/(n-1))*I_right - ((BBB.r*big.size())/(n-1))*I_right - ((BB.r*big.size())/(n-1))*I_right - ((B.r*big.size())/(n-1))*I_right - ((C.r*big.size())/(n-1))*I_right;
            
            if(gain > bestgain){
                bestgain = gain;
                ret.bestsplitval = splitval;
                ret.bestAttr = dataTemp[0][i];
            }
        }
    }
    return ret;
}
node* DTL(string data[][6],int minleaf,int n,int maxleaf){
     node *root = new node(false,"",0,"");
    // for(int i = 0; i < n; i ++){
    //     for(int j = 0; j < 6; j ++){
    //       cout << data[i][j] << "   ";
    //     }
    //     cout << '\n';
    // }
    //  cout << '\n';
    dt rec;
    if(is_Null(data)){
        root->isLeaf = true;
        root->label = "Null";
        return root;
    }
    if(is_leaf(data,n)){
        root->isLeaf = true;
        root->label = data[1][5];
        return root;
    }
    if(n < maxleaf+1){
        string tmp = cal_majority(data,n);
        root->isLeaf = true;
        root->label = tmp;
        return root;
    }
    rec = choose_split(data,n); 
    int pos = attr_pos(rec.bestAttr);
    root->isLeaf = false;
    root->attr = rec.bestAttr;
    root->splitval = rec.bestsplitval ;
    vector<string> left;
    vector<string> right;
    for(int i = 0; i < n-1; i ++){
        if(stod(data[i+1][pos]) <= rec.bestsplitval){
            left.push_back(data[i+1][pos]);
        }else{
            right.push_back(data[i+1][pos]);
        }
    }
    string rightData[right.size()+1][6];
    string leftData[left.size()+1][6];
    int lc = 1;
    int rc = 1;
    for(int j =0; j< n; j++){
        if(j == 0){
            for(int i = 0; i < 6; i++){
                leftData[j][i] = data[j][i];
                rightData[j][i] = data[j][i];
            }
        }else{
            if(stod(data[j][pos]) <= rec.bestsplitval){
                for(int i = 0; i < 6;i ++){
                    leftData[lc][i] = data[j][i];
                }
                lc++;
            }else{
                for(int i = 0; i < 6;i ++){
                    rightData[rc][i] = data[j][i];
                }
                rc++;
            }
        }
    }
     if (minleaf < maxleaf){
            root->left = DTL(leftData,minleaf+1,left.size()+1,maxleaf);
            root->right = DTL(rightData,minleaf+1,right.size()+1,maxleaf);
     }
     return root;
}
void predict(string data[], node* rootT){
    while(rootT->isLeaf == false){
        int t = attr_pos(rootT->attr);
        if( stod(data[t]) <= rootT->splitval){
            //cout << data[t] << " <= " << rootT->splitval << endl;
            rootT = rootT->left;
        }else{
            //cout << data[t] << " > " << rootT->splitval << endl;
            rootT = rootT->right;
        }
    }
    cout << rootT->label;
    cout << endl;
}
int main(int argc, char *argv[])
{   
    node *root;
    std::ifstream infile(argv[1]);
    std::ifstream file(argv[1]);

    std::ifstream testC(argv[2]);
    std::ifstream test(argv[2]);
    int minleaf = 1;
    std::string line;
    std::string tmp;
    std::string lineT;
    std::string tmpT;
    bool isAttr = 1;
    int iterator = 0;
    int counter = 0;
    std::string mt = argv[3];
    int maxleaf = stoi(mt);
    while (std::getline(infile, tmp)){counter ++;}
    string data[counter][6];
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        string attribute[6] = {"WC_TA", "RE_TA", "EBIT_TA", "MVE_BVTD", "S_TA", "Rating"};
        if (!(iss >> attribute[0] >> attribute[1] >> attribute[2] >> attribute[3] >> attribute[4] >> attribute[5])){
             break; 
        }
        if(isAttr == 1){
            for(int i = 0; i < 6; i ++){
                data[iterator][i] = attribute[i];
            }
            isAttr = 0;
            iterator++;
        }else{
            for(int i = 0; i < 6; i ++){
                data[iterator][i] = attribute[i];
            }
            iterator ++;
        }
    }
    root = DTL(data,minleaf,counter,maxleaf);
    int ct = 0;
    iterator = 0;
    isAttr = 1;
    while (std::getline(testC, tmpT)){ct ++;}
    string dataT[ct][5];
    while (std::getline(test, lineT))
    {
        std::istringstream iss(lineT);
        string attribute[5] = {"WC_TA", "RE_TA", "EBIT_TA", "MVE_BVTD", "S_TA"};
        if (!(iss >> attribute[0] >> attribute[1] >> attribute[2] >> attribute[3] >> attribute[4])){
             break; 
        }
        if(isAttr == 1){
            for(int i = 0; i < 5; i ++){
                dataT[iterator][i] = attribute[i];
            }
            isAttr = 0;
            iterator++;
        }else{
            for(int i = 0; i < 5; i ++){
                dataT[iterator][i] = attribute[i];
            }
            iterator ++;
        }
    }
    for(int i = 1; i < ct ; i ++){
        string temp[5];
        for(int j = 0; j < 5; j ++){
            temp[j] = dataT[i][j];
        }
        predict(temp,root);
    }
}