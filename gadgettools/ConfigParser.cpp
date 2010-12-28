#include <stdlib.h>

#include "ConfigParser.h"


bool ConfigParser::findSection(string section)
{
	return (findSectionPosition(section) != size_t(NOTFOUND))? true : false;
}

string ConfigParser::getStringVal(string section, string param)
{
	size_t start_sec = config.find("[" + section + "]");
	size_t hit = 0, hitend = 0;
	string value = "";
	
	if (start_sec == string::npos)
		return "";
	
	hit = config.find("\n" + param + "=", start_sec + 1);
	hitend = config.find("\n", hit + 1, 1);
	
	if (hit == string::npos || hitend == string::npos)
		return "";
	
	// Take into account brackets
	hit = hit + param.size() + 2;
	
	return config.substr(hit, hitend - hit);
}

size_t ConfigParser::getIntVal(string section, string param)
{
	return atoi(getStringVal(section, param).c_str());
}


float ConfigParser::getFloatVal(string section, string param)
{
	return atof(getStringVal(section, param).c_str());
}

void ConfigParser::add(string section, string param, string value)
{
	if (findSection(section))
	{
		// Section exists
		size_t insertLocation = findSectionPosition(section) + section.size() + 3;
		config.insert(insertLocation, param + "=" + value + "\n");
	}
	else
	{
		// Must add section
		config.append("\n[" + section + "]\n" + param + "=" + value + "\n");
	}
}

void ConfigParser::add(string section, string param, double value)
{
	stringstream out;
	out << value;
	string val = out.str();
	this->add(section, param, val);
}

void ConfigParser::add(string section, string param, size_t value)
{
	stringstream out;
	out << value;
	string val = out.str();
	this->add(section, param, val);
}

size_t ConfigParser::findSectionPosition(string section)
{
	size_t pos = 0;
	
	if (config.substr(0, section.size() + 3).compare("[" + section + "]\n") == 0)
		return pos;
	else
	{
		pos = config.find("\n[" + section + "]\n");
		
		if (pos == string::npos)
			return NOTFOUND;
		else
			return pos + 1;
	}	
}
