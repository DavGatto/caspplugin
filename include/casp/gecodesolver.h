#ifndef GECODESOLVER_H_
#define GECODESOLVER_H_

#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <string>
#include <map>

#include <boost/shared_ptr.hpp>

#include "simpleparser.h"


class GecodeSolver: public Gecode::MaximizeSpace {
public:
	/// Actual model
	GecodeSolver(std::vector<std::string> sumData, std::string domain, std::string globalConstraintName, std::string globalConstraintValue,
			boost::shared_ptr<SimpleParser> simpleParser);

	GecodeSolver(bool share, GecodeSolver& s) :
			Gecode::MaximizeSpace(share, s),
			_minValue(s._minValue),
			_maxValue(s._maxValue),
			_constraintVariables(s._constraintVariables),
			_costVariable(s._costVariable),
			_simpleParser(s._simpleParser) {
	}

	void propagate(std::vector<std::string> expressions);
	void propagate(std::string expression);

	virtual Gecode::Space*
	copy(bool share) {
		return new GecodeSolver(share, *this);
	}
	virtual Gecode::IntVar cost() const {
		return _costVariable;
	}
private:
	int _minValue;
	int _maxValue;

	std::vector<std::string> _sumData;

	boost::shared_ptr<SimpleParser> _simpleParser;

	std::map<std::string, Gecode::IntVar> _constraintVariables;
	Gecode::IntVar _costVariable;

	Gecode::LinExpr makeExpression(ParseTree* tree);
};



#endif /* GECODESOLVER_H_ */
