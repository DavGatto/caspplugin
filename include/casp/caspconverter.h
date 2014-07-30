#ifndef CASPCONVERTER_H_
#define CASPCONVERTER_H_

#include "dlvhex2/PluginInterface.h"
#include <utility>

class DefaultConverter : public dlvhex::PluginConverter {
public:
	DefaultConverter();
	virtual ~DefaultConverter();

	virtual void convert(std::istream& i, std::ostream& o);
};

bool syntaxCheckCnsAtm(std::string s, int LineNum);

class CaspConverter : public dlvhex::PluginConverter {
public:
	CaspConverter();
	virtual ~CaspConverter();

	virtual void convert(std::istream& i, std::ostream& o);

};



#endif /* CASPCONVERTER_H_ */
