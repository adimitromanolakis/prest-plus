/*
 * options.cpp
 *
 *  Created on: Nov 3, 2008
 *      Author: Apostolos Dimitromanolakis
 */

#include "options.h"

using namespace std;

options::options(int _argc, char **_argv)
{
    for (int i = 0; i < _argc; i++)
	argv.push_back(_argv[i]);

    argc = argv.size();

    used.resize(argc + 1);
    for (int i = 0; i < argc + 1; i++)
	used[i] = 0;

    // Read in script file if option present
    if (get("--script", 1)) {
	readscript(arg(0));
    }

}

void options::readscript(char *filename)
{
    ifstream in;
    in.open(filename);

    if (!in) {
	cout << "Error: Cannot open script file '" << filename << "'.\n";
	exit(1);
    }

    string s;

    while (in >> s) {
	argv.push_back(s);
	used.push_back(0);
    }

    argc = argv.size();
}

/*
 * Check for option opt_s, taking num_args arguments and return 1 if found
 * Also store arguments to arg array.
 */
int options::get(string opt_s, int num_args)
{
    int found = 0;

    const char *opt = opt_s.c_str();

    // find the option in the command line
    for (int i = 1; i < argc; i++) {
	if (!used[i]) if (strcmp(argv[i].c_str(), opt) == 0) {
	    found = 1;

	    used[i] = 1;

	    if (i + num_args >= argc) {
		cout << "Error: Option " << opt << " requires " << num_args << " arguments.\n";
		exit(1);
	    }

	    // copy the arguments to the arguments[] array
	    int k = 0;
	    i++;

	    while (k < num_args) {
		used[i] = 1;
		arguments[k++] = (char*) argv[i++].c_str();
	    }

	    i--;
	}

    }

    return found;
}

/*
 * Print all the command line options
 */
void options::print()
{
    for (int i = 1; i < argc; i++) {
	if (argv[i] == "--script") {
	    i++;
	    continue;
	}
	cout << argv[i] << " ";
    }
}

/*
 *  Check for unused options and report errors.
 */
void options::cleanup()
{
    // check for unused options
    int ok = 1;

    for (int k = 1; k < argc; k++)
	if (!used[k]) {
	    cout << "Error: Option " << argv[k] << " not understood.\n";
	    ok = 0;
	    break;
	}

    if (ok == 0) {
	cout << "\nTry using --help for a list of command line options.\n";
	exit(1);
    }
}

options::~options()
{

}
