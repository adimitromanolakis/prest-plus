
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <map>
#include <sstream>
#include <string>

using namespace std;

string join(const vector<string> &f, string delim = " ")
{
	int i;
	string s = "";
	if(f.size() > 0) { for(i=0;i<f.size()-1;i++) { s += f[i]; s += " "; }
	s  += f[f.size()-1]; }
	return s;
}


vector<string> split(string s)
{
	istringstream in(s);
	vector <string> ret;
	string a;

	while(in >> a) ret.push_back(a);
	return ret;
}

vector<string> split(string s, char delim )
{
	istringstream in(s);
	vector <string> ret;
	string a;

	while(getline(in, a, delim) >> a) ret.push_back(a);
	return ret;
}

 
template<class X>
void add(const vector<X> &a, const vector<X> &b)
{

}

template<class X>
vector<X> sub(vector<X> &a, int l, int u) {
   if(u==-1) u = a.size()-1;
    vector <X> t;
    if(l < a.size()-1 && u>=l)
    for(int i=l;i<=u;i++) t.push_back(a[i]);
    return t;
}


string conv(double x) { char a[20]; sprintf(a,"%f",x); return (string)a; }
string conv_int(int x) { char a[20]; sprintf(a,"%d",x); return (string)a; }

