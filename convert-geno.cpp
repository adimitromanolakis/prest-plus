
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include <sstream>
#include <map>
#include <algorithm>

#include <math.h>
#include "utils.cpp"

using namespace std;

class person {
public:
    string fid,id,dad,mom,sex,pheno;

    string geno;
    int nochrom;
    person() { nochrom = 0; }
};


vector<person> ped;
map<string,int> pindex;


vector<string> nextline(ifstream &i) {
    string s;
    getline(i,s);
    return split(s);
}

char str1[10]; 

vector< vector<string> > mapdata;


double mapthetatocm(double theta, int maptype)
{
  double cm = 0;
 
  /* no-interference model */
  if ( maptype == 1 ) {
    cm = -0.5*log(1-2*theta)*100;
  }
  /* kosambi model */
  else if ( maptype == 2 ) {
    cm = 0.25*log((1+2*theta)/(1-2*theta))*100;
  }

  return cm;
}

main(int argc,char *argv[])
{
    if(argc!=3) {
        cout << "convert-geno: A tool to convert old linkage format files to prest-plus file format\n\n";
        cout << "Usage: convert-geno pedigree_file index_file \n";
        cout << " where index file contains one line per chromosome with map and ped files:\n"
                "  map_file ped_file\n\n";

        exit(0);
    }

    cout << "Converting old format files to new format...\n\n";
    
    string pedigree_file = argv[1];
    string chrom_files = argv[2];

    ifstream in (pedigree_file.c_str());
    ifstream chrom (chrom_files.c_str());

    cout << "Reading pedigree file: " << argv[1] << "\n";


    person p;
    string s;
    
    while(getline(in,s)) {
        vector<string> v = split(s);
        p.id = v[1];
        p.fid = v[0];
        p.dad = v[2];
        p.mom = v[3];
        p.sex = v[4];
        p.pheno = v[5];
        string ind = p.fid + ".." + p.id;
        pindex[ind] = ped.size();
        ped.push_back(p);
        sub(v,0,1);
        
        // cout << join(sub(v,0,3)) << "\n";
    }

    int CHR = 0;
    string fped_s,fmap_s;
    
    while(chrom >> fmap_s >> fped_s) {
        vector<string> v;
        CHR ++;
        cout << "Reading chromosome files: " << fmap_s << " " << fped_s << "\n";
        ifstream fmap(fmap_s.c_str());
        ifstream fped(fped_s.c_str());

        if(fmap_s==(string)"-") cout << "no map file!" << "\n";

        int nloci, t1,t2,t3;

        v = nextline(fped);
        nloci = (v.size() - 6)/2;
        cout << "    #loci = " << nloci << "\n";

        fped.seekg(0);
        
        while(getline(fped,s)) {
            v = split(s);
            string ind = v[0] + ".." + v[1];

            double N = (v.size() - 6)/2;
            if(N != nloci) cout << "Error: individual %s %s in ped file %s has inconsistent # of markers"
                            << N << " instead of " << nloci << "\n";
            
            if(pindex.find(ind) != pindex.end()) {
                int k = pindex[ind] ;

                if(ped[k].nochrom == CHR) {
                    cout << "Error in pedigree file: duplicate individual " << ind << "\n";
                    cout << " -> skipping\n";
                    continue;

                }
                //string &s = ped [ ].geno;
                ped[k].geno  += "      " + join(sub(v,6,-1));
                ped[k].nochrom = CHR;

                //cout << ind << " " <<  k << " " << v.size() <<  " " << ped[k].geno.size() << "  MM\n";

            }
        }

        for(int i=0;i<ped.size();i++) if(ped[i].nochrom != CHR) {
            cout << "Warning: filling zeroes for missing individual: " << ped[i].fid << " " << ped[i].id << "\n";

            vector<string> t(2*nloci);
            for(int z=0;z<2*nloci;z++) t[z] = "0";
            ped[i].geno += ("      " + join(t));
        }

        int havemap = (fmap_s != "-");

        if(havemap) {
            v = nextline(fmap);
            nloci = atoi(v[0].c_str());
        
            //cout << "map file loci = " << nloci << "\n";
        }

        int pos_start = mapdata.size();
        
        for(int i=0;i<nloci;i++) {
            char name[100];
            sprintf(name,"CHR%d.M%d", CHR, i+1);

            
            vector<string> mapl(4);
            mapl[0] = conv_int(CHR);
            mapl[1] = name;
            mapl[2] = "0";
            mapl[3] = "0";

            if(havemap) {
                v = nextline(fmap);
                int noalle = atoi(v[1].c_str());

                v = nextline(fmap);

                mapl.push_back(conv_int(noalle));
                for(int j=0;j<noalle;j++) mapl.push_back(v[j]);
            }

            mapdata.push_back(mapl);
        }

        if(havemap) {
            vector<double> recomb;
            v = nextline(fmap);

            double cm = 0;

            mapdata[pos_start][2] = conv(cm);

            for(int i=0;i<nloci-1;i++) {
                cm += mapthetatocm(atof(v[i].c_str()), 2);
                mapdata[pos_start+i+1][2] = conv(cm);
            }
        }

        

    }
   

    ofstream o("prest.map");

    for(int i=0;i<mapdata.size();i++)
        o << join( mapdata[i] ) << "\n";


    ofstream fped("prest.ped");



    for(int i=0;i<ped.size();i++)
    {
        fped << ped[i].fid << " " << ped[i].id << " " 
             << ped[i].dad << " " << ped[i].mom << " " 
             << ped[i].sex << " " << ped[i].pheno << " "
             << ped[i].geno << "\n";

        // cout << ped[i].geno.size() << " " <<  split(ped[i].geno).size() << "\n";

    }


    cout << "\n\n"
            "Finished conversion.\n"

            "Total individuals: " << ped.size() << "\n"
            "Total markers: " << mapdata.size() << "\n\n"
            "Output files are prest.map and prest.ped\n\n"
    "You can run prest-plus now by typing:\n  prest --file prest.ped --map prest.map --wped --aped\n\n"
            ;

}
