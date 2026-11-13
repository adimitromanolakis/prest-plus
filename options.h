/*
 * options.h
 *
 *  Created on: Nov 3, 2008
 *      Author: vaio
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

class options
{
private:

	int argc;
	vector<string> argv;
	vector<int> used;

	char *arguments[10];

public:
	options(int _argc, char **_argv);
	virtual ~options();

	void cleanup();
	int get(string opt, int num_args);
	void readscript(char *filename);
	void print();

	char * arg(int i) { return arguments[i]; }
};

#endif /* OPTIONS_H_ */
