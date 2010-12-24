#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <string>
#include <iostream>
#include <sstream>

#define NOTFOUND -1

using namespace std;

class ConfigParser 
{	
public:
	ConfigParser() { config = ""; }
	
	void parse(string inc) { config = inc; }
	void parse(char * inc) { config = inc; }
	
	bool findSection(string section);
	
	string getStringVal(string section, string param);
	size_t getIntVal(string section, string param);
	float getFloatVal(string section, string param);
	
	void add(string section, string param, string value);
	void add(string section, string param, double value);
	void add(string section, string param, size_t value);
	
	string exportConf() { return config; };
	const char * exportCharArray() { return config.c_str(); }
	
private:
	string config;

	size_t findSectionPosition(string section);
};

#endif //CONFIGPARSER_H
