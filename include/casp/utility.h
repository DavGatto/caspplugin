#ifndef UTILITY_H_
#define UTILITY_H_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <set>


#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <dlvhex2/Nogood.h>
#include <dlvhex2/Term.h>

/**
 * @brief This helper method is used to check whether string is a variable
 */
bool isVariable(std::string s);

bool isNumber(std::string s);

int numofdgts(int n);

/**
 * @brief This helper method removes quotes from the string
 */
std::string removeQuotes(std::string s);

/**
 * @brief This method is used to change the operator to invertible one
 * in "not_expr" terms, e.g. not_expr("3>5") becomes expr("3<=5").
 */
std::string replaceInvertibleOperator(std::string expr);

/**
 * @brief This helper method checks, whether string is a separator
 */
bool isSeparator(std::string s);

/**
 * @brief This helper method get's value of string, which is
 * contained inside of the brackets
 */
std::string getValueInsideBrackets(std::string s);

/** Used in caspconverter: */

/**
 * @brief This method takes a predicate constraint variable and outputs
 * the ASP variables which appear in it as a comma separated string; also
 * it increments the int arity of the number of ASP variables it finds
 */
std::string extractVariables(std::string s, int &arity);

/**
 * @brief This method changes in '?' all the ASP variables in its input
 */
std::string generalizeVariables(std::string s);

/**
 * @brief This method returns true iff its input is + | - | * | /
 */
bool isBinOperator(char s);

/**
 * @brief This method returns true iff its input is = | ! | < | >
 */
bool isComparison(char s);
/**
 * @brief This method returns true iff its input is a valid symbol in constraint atom syntax
 */
bool isValidConstraintSym(char s);

/*int get_lower_bound(std::string domain);

int get_upper_bound(std::string domain);*/

/**
 * @brief This method returns < if input is > and viceversa, and returns its input unchanged otherwise
 */
char invertComparison(char c);

/**
 * @brief This method applies invertComparison to all the comparison symbols in its input
 */
std::string invertComparison(std::string s);

/** Canonizer Methods */

/**
 * @brief This enum lists the possible kinds of comparison that can appear in an expression
 */
enum ComparisonType {
	EQ, NEQ, GT, LT, GEQ, LEQ
};

/**
 * @brief This struct represents a linear expression
 */
struct expression{
    public:
        std::string ID; /** The unique name of the expression, e.g. not_expr_5_3("age(X)+Y=Z",X,Y,Z) */
        std::string raw = ""; /** The form in which the expression appears */
        std::set<int> scope; /** The canonizer-variables which appear in the expression */
        std::string lhs = ""; /** The LHS of the canonical form of the expression */
        std::string RHS = ""; /** The RHS as a string */
        bool intrhs = false; /** This is set to true if rhs can be expressed as an int */
        int rhs = 0; /** The RHS of the canonical form of the expression, if it is an int */
        ComparisonType co; /** the comparison type of the expression */
        bool hasmin = false; /** this is true if a min val is known */
        int minval = 0; /** the min val */
        bool hasmax = false; /** this is true if a max val is known */
        int maxval = 0; /** the max val */
};

/*struct signednogood{
    bool sgn = false;
    dlvhex::Nogood ng;
};*/

/**
 * @brief This method returns the abs difference between maxval and minval of its input, if they are known.
 * Otherwise cerr a message and outputs -1
 */
int chances(expression e);

/**
 * @brief This method changes minval and maxval of two expressions having the same lhs
 * It changes the min and max of the first one wrt those of the second one, so normally is used twice
 * ie relate(a,b); realte(b,a);
 */
void relate(expression a, expression b);

/**
 * @brief This struct describes fractions
 */
struct fraction{
    int numer=1;
    int denom=1;
};

/**
 * @brief This method simplifies a fraction
 */
fraction reduce(fraction n);

/**
 * @brief This method defines the sum of fractions
 */
fraction operator+(fraction a, fraction b);

/**
 * @brief This method returns the fraction wrote in the input string
 */
fraction atofraction(std::string s);

/**
 * @brief This method returns the number found immediately after position x
 * it is used with canonizer variables like X23 or V12
 */
int takevarnum(std::string t, int x);

/**
 * @brief This method defines an ordering for canonizer variables
 */
bool variableorder(std::string i, std::string j);

/**
 * @brief This method defines an order for compound terms like -2/3*X7*V6/X8
 */
bool termorder(std::string i, std::string j);

/**
 * @brief This method defines an order for constraint variables
 */
bool ConstrVarOrder(std::string i, std::string j);
/**
 * @brief This class is used to deal with constraints atoms rewriting them in a canonical form
 */

/**
* @brief This method returns true iff s is made only of digits and parenthesis
*/
bool isNumberOrPar(std::string s);

class Canonizer{


    /**
    * @brief This method takes an expression string s and returns a string in which constraint variables are replaced by canonizer variables
    * and stores in varIDs the (original var)<->(canonizer var) relations, and in scope the scope of the expression
    */
    std::string makeVariablesIDs(std::string s, std::set<int> &scope);

    /**
    * @brief This method
    */
    std::string aliasDenum(std::string s);

    /**
    * @brief This method
    */
    std::string simplifyTerm(std::string s);

    /**
    * @brief This method
    */
    std::string sortTerm(std::string s);

    /**
    * @brief This method
    */
    std::string computeAsString(std::vector<std::string> v);

    /**
    * @brief This method
    */
    std::string simplify(std::string s,std::string oper);

    /**
    * @brief This method
    */
    std::string sortExpr(std::string s);

    /**
    * @brief This method
    */
    std::string computeSums(std::string s);

    /**
    * @brief This method
    */
    std::string multiplyByMinusOne(std::string s);

    /**
    * @brief This method
    */
    std::string takeToRHS(std::string s);

    /**
    * @brief This method
    */
    std::string makeLHSonly(std::string s);

    /**
    * @brief This method
    */
    expression makeExpression(std::string in,std::string id);

public:
    Canonizer();

    std::map<std::string, int> varIDs; /** Stores the original names of the canonized variables */

    std::map<int,std::string> invIDs; /** Just the inverse map of the former */

    int freeASPvarnum = 0; /** Number of free ASP vars used as constraint vars */

    std::vector<expression> Expressions; /** stores canonized Expressions */

    int dom_min = 0, dom_max = 0;

    double prob_of_is = 0;

    /**
    * @brief This method returns the canonical form of the expr string in its input, and stores its scope in scope
    */
    std::string canonicalForm(std::string s, std::set<int> &scope);

    /**
    * @brief This method returns the canonical form of the expr string in its input, without storing the scope
    */
    std::string canonicalForm(std::string s);

    /**
    * @brief This method returns a string in which canonizer vars are replaced by their original names
    */
    std::string backToVarNames(std::string s);

    /**
    * @brief This method returns the generic form of a nonground constraint atom,
    * e.g.  expr_3_2("age(X)$=height(Y)",X,Y)  ->  expr_3_2("age(Var1)$=height(Var2)",Var1,Var2)
    */
    std::string generalizeNongroundExpr(std::string s);

    /**
    * @brief This method
    */
    std::string generalizeVars(std::string s, int arity);

    /**
    * @brief This method adds the expression described by string s in the vector<expression> Expressions
    * assigning to it the string id as its .ID
    */
    void addExpression(std::string s,std::string id);

    /**
    * @brief This method returns a DLVHEX rule if two expressions in Expressions are found to be incompatible
    * with each other, or one implied by the other and so on.
    */
    std::string analize();


    //std::string nonground_analize();
    /**
    * @brief This method couts the Expression[k]{.lhs, comparisontype, rhs} from the whole vector Expressions
    */
    void printExpressions();

    /**
    * @brief This method couts the varIDs content
    */
    void printvarIDs();

    /**
    * @brief This method clears varIDs, Expressions and freeASPvarnum
    */
    void clear();

};


#endif /* UTILITY_H_ */
