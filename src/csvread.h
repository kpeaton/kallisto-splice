#ifndef CSVREAD_H
#define CSVREAD_H
// Adapted from a Stack Overflow answer from user Loki Astari: http://stackoverflow.com/a/1120224/52738

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class CSVRow
{
public:
	std::string const& operator[](std::size_t index) const;
	std::size_t size() const;
	void readNextRow(std::istream& str);
private:
	std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data);

#endif // CSVREAD_H