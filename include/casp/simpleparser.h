#ifndef SIMPLEPARSER_H_
#define SIMPLEPARSER_H_

#include <string>
#include <map>
#include <set>


enum ElementType {
	NUMBER, CONSTRAINT_VARIABLE, OPERATOR
};

struct ParseTree {
	std::string value;
	struct ParseTree *left, *right;
	ElementType type;
};

const int MAX = 500;

class SimpleParser {
private:
	std::map<std::string, ParseTree> _cachedTrees;

	int isOperator(std::string s);
	int priority(std::string s);
	std::string convertToPostfix(std::string infix);
public:
	SimpleParser();
	std::set<std::string> getConstraintVariables(struct ParseTree* root);
	void makeTree(std::string infix, struct ParseTree** root);
	void deleteTree(struct ParseTree* root);
};

#endif /* SIMPLEPARSER_H_ */
