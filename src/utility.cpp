#include "casp/utility.h"

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

using namespace std;

/**
 * @brief This helper method is used to check whether string is a variable
 */
bool isVariable(string s) {
	bool res = true;
	if (s[0] >= 'A' && s[0] <= 'Z') {
		for (int i = 1; i < s.length(); i++) {
			if (!isalnum(s[i]) && s[i] != '_') {
				res = false;
				break;
			}
		}
	} else
		res = false;
	return res;
}

bool isNumber(string s)
{
    for(int i=0;i<s.length();i++)
        if(!isdigit(s[i]))
            return false;
    return true;
}

int numofdgts(int n)
{
    int p=10, d=1;
    while(true){
        if(abs(n)<p)
            return d;
        else{
            d++; p *= 10;
        }
    }
}

string removeQuotes(string s) {
	return s.substr(1, s.length() - 2);
}
/**
 * This method is used to change the operator to invertible one
 * in "not_expr" terms, e.g. not_expr("3>5") becomes expr("3<=5").
 */
string replaceInvertibleOperator(string expr) {
	typedef pair<string, string> string_pair;

	vector<string_pair> invertibleOperators;

	invertibleOperators.push_back(make_pair("==", "!="));
	invertibleOperators.push_back(make_pair("!=", "=="));
	invertibleOperators.push_back(make_pair(">=", "<"));
	invertibleOperators.push_back(make_pair("<=", ">"));
	invertibleOperators.push_back(make_pair(">", "<="));
	invertibleOperators.push_back(make_pair("<", ">="));

	BOOST_FOREACH (string_pair invertibleOperator, invertibleOperators) {
		if (expr.find (invertibleOperator.first) != string::npos) {
			boost::replace_all(expr, invertibleOperator.first, invertibleOperator.second);
			break;
		}
	}

	return expr;
}

bool isSeparator(string s)
{
	return s == "," || s == "v" || s == "." || s == ":" || s == "not";
}

string getValueInsideBrackets(string s) {
	int startIndex = s.find("(") + 1;
	int endIndex = s.find(")");
	return s.substr(startIndex, endIndex - startIndex);
}

/** Used in caspconverter.cpp : */

string extractVariables(string s, int &arity)
{
    string v = "";
    set<string> var;
    if(isVariable(s)){arity++; return s;}
    int b = 0, l = 0;
    while(b<s.length()){
        while( !isupper(s[b]) ){ b++; if(b == s.length()){ break;} }
        if(b>=s.length()){ break; }
        if(v != ""){ v.append(","); }
        l = 0;
        while( isalnum(s[b+l]) || s[b+l] == '_' ){ if(b+l >= s.length()){ break;} l++; }
        if( isVariable(s.substr(b,l)) ){ if(var.count(s.substr(b,l))==0){v.append(s.substr(b,l)); arity++;} }
        if(b+l >= s.length()) { break;}
        b+=l;
    }
    return v;
}

string generalizeVariables(string s)
{
    bool var = 0, v=0;
    string o="";
    for(int i=0;i<s.length();i++){
        if(( i == 0 || !isupper(s[i-1]) ) && isupper(s[i])){ var = 1;}
        if(!isalnum(s[i]) && s[i]!='_' && !isBinOperator(s[i])){ var = 0;}
        if(var){ s[i] = '?';}
        if(s[i]=='?'){ if(i==0 || s[i-1]!='?') { o.append("?");}}
        else { o.append(s.substr(i,1));}
    }
    return o;
}

bool isBinOperator(char s)
{
    return s == '+' || s == '-' || s == '*' || s == '/';
}

bool isComparison(char s)
{
    return s == '=' || s == '!' || s == '<' || s == '>';
}

bool isValidConstraintSym(char s)
{
    if(isalnum(s)) { return true; }
    return s == '=' || s == '!' || s == '<' || s == '>' || s == '+' || s == '-' || s == '*' || s == '/' || s == '$' || s == '}' || s == '{' || s == ';' || s == '_' || s=='\0';
}

/*int get_lower_bound(string domain)
{
    int i=0;
    while(isdigit(domain[i]))
        i++;
    return atoi(domain.substr(0,i).c_str());
}

int get_upper_bound(string domain)
{
    int i=domain.length()-1;
    while(isdigit(domain[i]))
        i--;
    return atoi(domain.substr(i+1).c_str());
}*/

/** Canonizer Methods */
string Canonizer::makeLHSonly(string s)
{
    //cout << "Called 'makeLHSonly' on '" << s << "'" << endl;
    stringstream op;
    op.str(string());
    int parlevel = 0;
    for(int i=s.length()-1;i>0;i--){
        if(isComparison(s[i])){
            if(s[i+1]!='+' && s[i+1]!='-'){
                if(isComparison(s[i-1])){
                    op << s[i-1] << s[i];
                    s[i] = '-';
                    s.erase(i-1,1);
                }
                else{
                    op << s[i];
                    s[i] = '-';
                }
            }
            else{
                if(isComparison(s[i-1])){
                    op << s[i-1] << s[i];
                    s.erase(i-1,2);
                }
                else{
                    op << s[i];
                    s.erase(i,1);
                }
            }
            s.append(op.str());
            s.append("0");
            break;
        }
        if(s[i]=='{') /** This is used in caspconverter after the switch ( -> { */
            parlevel++;
        if(s[i]=='}')
            parlevel--;
        if(parlevel==0){
            if(s[i]=='+')
                s[i]='-';
            else if(s[i]=='-')
                s[i]='+';
        }
    }

    return s;
}

char invertComparison(char c)
{
    if(c=='>')
        return '<';
    else if(c=='<')
        return '>';
    else
        return c;
}

string invertComparison(string s)
{
    for(int i=0;i<s.length();i++)
        if(isComparison(s[i]))
            s[i] = invertComparison(s[i]);
    return s;
}

using namespace boost;

fraction reduce(fraction n)
{
    fraction r;
    if(n.numer*n.denom==-1 && n.denom<0){
        n.numer *= -1;
        n.denom *= -1;
    }
    if(n.numer%n.denom==0){
        r.numer = n.numer / n.denom;
        r.denom = 1;
        return r;
    }
    else{
        int i = std::min(n.numer,n.denom)-1;
        while(i>1){
            if(n.numer%i==0 && n.denom%i==0){
                n.numer /= i;
                n.denom /= i;
            }
            i--;
        }
    }
    return n;
}

fraction operator+(fraction a, fraction b)
{
    fraction r;
    r.denom = a.denom * b.denom;
    r.numer = b.denom*a.numer + a.denom*b.numer;
    return reduce(r);
}

fraction atofraction(string s)
{
    fraction r;
    size_t found = s.find('/');
    if(found != string::npos){
        r.numer = atoi(s.substr(0,found).c_str());
        r.denom = atoi(s.substr(found+1).c_str());
    }
    else{
        r.numer = atoi(s.c_str());
    }
    if(r.denom==0){
        cerr << "atofraction problem!" << endl;
        r.denom = 1;
        return r;
    }
    return reduce(r);
}

void relate(expression a, expression b)
{
    if(a.lhs!=b.lhs)
        return;

    if(!b.hasmin && !b.hasmax)
        return;

    if(b.hasmin){
        if(a.hasmin){
            if(a.minval<b.minval)
                a.minval = b.minval;
        }
        else{
            a.hasmin = true;
            a.minval = b.minval;
        }
    }
    if(b.hasmax){
        if(a.hasmax){
            if(a.maxval>b.maxval)
                a.maxval = b.maxval;
        }
        else{
            a.hasmax = true;
            a.maxval = b.maxval;
        }
    }
}

int chances(expression e)
{
    if(e.hasmin && e.hasmax)
        return abs(e.maxval - e.minval);
    else{
        cerr << "Called 'chances' on unbounded expression." << endl;
        return -1;
    }

}

int takevarnum(string t, int x)
{
    int l = 1;
    while(isdigit(t[x+l]))
        l++;
    return atoi(t.substr(x+1,x+l-1).c_str());
}

bool variableorder(string i, string j)
{

    size_t xi = i.find("X"), xj = j.find("X"),
           vi = i.find("V"), vj = j.find("V");
    if(xi==string::npos){ /** If i has no Xs... */
        if(xj==string::npos){ /** and j has no Xs... */
            if(vi==string::npos){ /** and i has no Vs... */
                if (vj==string::npos){ /** and j has no Vs... */
                    if(i[0]=='/')
                        return false; /** then i>=j if i begins with / */
                    else
                        return true; /** otherwise i<j */
                }
                else
                    return true; /** i<j as j has Vs and i has not */
            }
            else{
                if(vj==string::npos)
                    return false; /** i>=j as i has Vs and j has not */
                else
                    return takevarnum(i,vi) < takevarnum(j,vj); /** i<j iff varnumber_i<varnumber_j */
            }
        }
        else
            return true; /** i<j as j has Xs and i has not */
    }
    else{
        if(xj==string::npos)
            return false; /** i>=j as i has Xs and j has not */
        else
            return takevarnum(i,xi) < takevarnum(j,xj); /** i<j iff varnumber_i<varnumber_j */
    }
}

bool termorder(string i, string j)
{
    /** Store first chars (which are + or -) in opi and opj, and erase from original strigs */
    int ix = 0, jx = 0, iv = 0, jv = 0;
    char opi = i[0], opj = j[0];
    i.erase(0,1); j.erase(0,1);
    /** Tokenize both i and j wrt * and /, dropping operations */
    char_separator<char> sep("*/ ", "", drop_empty_tokens);

	tokenizer<char_separator<char> > itokens(i, sep);
	tokenizer<char_separator<char> > jtokens(j, sep);

	vector<string> ivars, jvars;

	for ( tokenizer<char_separator<char> >::iterator it = itokens.begin(); it != itokens.end(); ++it) {
        string tok = *it;

        if(tok[0]=='X'){
            ix++;
            ivars.push_back(tok);
        }
        if(tok[0]=='V'){
            iv++;
            ivars.push_back(tok);
        }
	}

	for ( tokenizer<char_separator<char> >::iterator it = jtokens.begin(); it != jtokens.end(); ++it) {
        string tok = *it;

        if(tok[0]=='X'){
            jx++;
            jvars.push_back(tok);
        }
        if(tok[0]=='V'){
            jv++;
            jvars.push_back(tok);
        }
	}

	/** ix = # of Xvars in i, vx = # of Vvars in i, same for j */
	if(ix+iv < jx+jv)
        return true; /** i<j if i has less variables */
    else if(ix+iv == jx+jv && ix < jx)
        return true; /** i<j if they have the same number of variables, but i has less Xs */

	/** Now, if i and j have the same number of Xs and Vs, then i<j if the first different variable between them is < according to variableorder */
	for(int h=0; h<ivars.size(); h++){
        if(variableorder(ivars[h],jvars[h]))
           return true;
        if(variableorder(jvars[h],ivars[h]))
           return false;
	}
	/** Now, if all the variables are the same between i and j, than we look at opi and opj */
    if(opi == '+' && opj == '-' )
        return true; /** '+' < '-' */
    else
        return false; /** in any other case, i >= j. */

}

bool ConstrVarOrder(string i, string j)
{
    if(isVariable(i)){
        if(isVariable(j))
            return i<j; /** i and j are both ASP var, return i<j in lexicographical sense */
        else
            return true; /** i is an ASP var while j is not, so i<j */
    }
    int ari = 0, arj = 0;
    string exi = extractVariables(i,ari), exj = extractVariables(j,arj);
    if(ari<arj)
        return true; /** i has smaller arity than j, so i<j */
    if(arj<ari)
        return false; /** j has smaller arity than i, so i>=j */
    size_t ip = i.find('('), jp = j.find('(');
    if(i.substr(0,ip+1)<j.substr(0,jp+1))
        return true; /** predicate name of i is < than pred name of j in lexicographical sense, so i<j */
    if(j.substr(0,jp+1)<i.substr(0,ip+1))
        return false; /** predicate name of j is < than pred name of i in lexicographical sense, so i>=j */
    /** Same predicate name, same arity. Let's check ASP vars names, beginning from the first, which is
        simply looking for the lexicographical ordering between exi and exj */
    return exi < exj;

}

bool isNumberOrPar(string s)
{
    for(int i=0;i<s.length();i++)
        if(!isdigit(s[i]) && s[i]!='(' && s[i]!=')')
            return false;
    return true;
}


Canonizer::Canonizer()
{
    varIDs.clear();
    Expressions.clear();
}



string Canonizer::makeVariablesIDs(string s, set<int> &scope)
{
    boost::replace_all(s,"{","(");
    boost::replace_all(s,"}",")");
    boost::replace_all(s,";",",");

    //cout << "Called 'makeVariablesIDs' '" << s << "'" << endl;
    scope.clear();
    //stringstream out;
    vector<string> variables;
    int i = 0;
    while(i<s.length()){
        if(isalpha(s[i])){
            int j=1, parlevel=0;
            while(isalnum(s[i+j]) || s[i+j]=='_' ){
                j++;
            }
            if(s[i+j]=='('){
                parlevel = 1;
                while(parlevel!=0){
                    j++;
                    if(s[i+j]=='(')
                        parlevel++;
                    if(s[i+j]==')')
                        parlevel--;
                }
                j++;

            }
            string var = s.substr(i,j+1);
            if(isBinOperator(var.back()))
                var.pop_back();
            variables.push_back(var);
            i += j;
        }
        i++;
    }

    std::sort(variables.begin(),variables.end(),ConstrVarOrder);
    int n = varIDs.size();
    for(int k=0;k<variables.size();k++){
        if(varIDs.count(variables[k])==0){
            varIDs.insert({variables[k],n}); n++;
        }
        size_t fv = s.find(variables[k]);
        while(fv!=string::npos){
            stringstream vvv;
            vvv << "X" << varIDs.at(variables[k]);
            s = s.substr(0,fv) + vvv.str() + s.substr(fv+variables[k].length());
            fv = s.find(variables[k]);

        }
    }

    for(auto it=varIDs.begin(); it!=varIDs.end(); it++){
            invIDs.insert({it->second,it->first});
    }

    return s;
}


string Canonizer::aliasDenum(string s)
{
    for(int i=0;i<s.length()-1;i++){
        if(s[i]=='/'&&s[i+1]=='X'){
            s[i] = '*';
            s[i+1] = 'V';
        }
    }
    return s;
}



string Canonizer::simplifyTerm(string s)
{
    for(int i=0;i<s.length();i++){
         if(s[i]=='X'||s[i]=='V'){
            stringstream var;
            if(s[i]=='X')
                var << "V" << takevarnum(s,i);
            else
                var << "X" << takevarnum(s,i);
            size_t found = s.find(var.str());
            if(found!=string::npos){
                int l = var.str().length();
                s.erase(found-1,l+1);

                if(i>0){
                    s.erase(i-1,l+1);

                }
                else{
                    s.erase(i,l+1);

                }
                if(i<s.length()){
                    return s.substr(0,i) + simplifyTerm(s.substr(i));
                }
                else
                    return s;
            }
         }
    }
    return s;
}


string Canonizer::sortTerm(string s)
{
    s = "*" + s;

    stringstream out;
    vector<string> tokenList;
    vector<int> vars;
    char_separator<char> sep(" ", "*/", drop_empty_tokens);
	tokenizer<char_separator<char> > tokens(s, sep);
	vector<int> tokenDeg;
    string op = "";
	for ( tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
        if(*it=="*"||*it=="/")
            op = *it;
		else{
            string token = op + *it;
            tokenList.push_back(token);
		}
	}

	sort(tokenList.begin(),tokenList.end(),variableorder);

	for (std::vector<string>::iterator it=tokenList.begin(); it!=tokenList.end(); ++it)
        out << *it;
    if(out.str()[0]=='*')
        return simplifyTerm(out.str().substr(1));
    if(out.str()[0]=='/')
        return simplifyTerm("1"+out.str());

    return simplifyTerm(out.str());
}

string Canonizer::computeAsString(vector<string> v)
{
    int n = 0;
    for(int i=0;i<v.size();i++){
        if(v[i][0]=='+')
            n += atoi(v[i].substr(1).c_str());
        else if(v[i][0]=='-')
            n -= atoi(v[i].substr(1).c_str());
    }
    stringstream o;
    o << n;
    return o.str();
}

string Canonizer::simplify(string s,string oper)
{

    s = "*"+s;
    int numer = 1, denom = 1;
    char_separator<char> sep(" ", "*/", drop_empty_tokens);
    tokenizer<char_separator<char> > tokens(s, sep);
	vector<string> tokenv;
    string op = "";
	for ( tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
        if(*it=="*"||*it=="/")
            op = *it;
        else{
            string N = *it;
            if(op=="*")
                numer *= atoi(N.c_str());
            else if(op=="/")
                denom *= atoi(N.c_str());
            op = "";
        }
	}

    int coef = numer/denom, remainder = numer - coef*denom;
    stringstream o;
    o << coef;
    if(remainder!=0)
        o << oper << remainder;

    return o.str();
}

string Canonizer::sortExpr(string s)
{
    if(s[0]!='-')
        s = "+" + aliasDenum(s);
    s = aliasDenum(s);
    string out = "";
    char_separator<char> sep(" ", "+-", drop_empty_tokens);
	tokenizer<char_separator<char> > tokens(s, sep);
	vector<string> tokenList;
    vector<int> tokenDeg;
    string op = "";
	for ( tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
		if(*it=="+" || *it=="-"){
            op = *it;
		}
		else{
            string token = *it;


            int deg = 0;
            for(int j=0;j<token.length();j++){
                if(token[j]=='X'||token[j]=='V'){
                    deg++;
                }
            }
            if(deg==0){
                token = simplify(token,op);
            }
            tokenDeg.push_back(deg);
            tokenList.push_back(op + sortTerm(token));
            op = "";

		}
	}
    auto it = max_element(std::begin(tokenDeg), std::end(tokenDeg));
    vector<string> samedeg;

    for(int d=0;d<=*it;d++){
        samedeg.clear();
        for(int i=0;i<tokenList.size();i++)
            if(tokenDeg[i]==d)
                samedeg.push_back(sortTerm(tokenList[i]));

        if(d==0){

            out.append(computeAsString(samedeg));
        }

        else{
            sort(samedeg.begin(),samedeg.end(),termorder);

            for (std::vector<string>::iterator it=samedeg.begin(); it!=samedeg.end(); ++it)
                out.append(*it);
        }
    }
    return out;
}

string Canonizer::computeSums(string s)
{
    s.append("+T");
    char_separator<char> sep(" ", "+-", drop_empty_tokens);
	tokenizer<char_separator<char> > tokens(s, sep);
	vector<string> tokenList;
	string previous = "", op = "";
	stringstream out;
	fraction prevcoef;
	for ( tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
        string token = *it;
        if(token=="+"||token=="-")
            op = token;
        else{
            int v = -1;
            for(int i=0;i<token.length();i++){
                if(token[i]=='X'||token[i]=='V'||token[i]=='T'){
                    v = i;
                    break;
                }
            }
            if(v==-1){
                out << op << token;
            }
            else{
                string term;
                fraction coef;
                if(v==0){
                    if(op == "-")
                        coef.numer *= -1;
                    term = token;
                }
                else{
                    string strcoef = op+token.substr(0,v-1);
                    coef = atofraction(strcoef);
                    term = token.substr(v);
                }
                if(previous==""){
                    previous = term;
                    prevcoef = coef;
                }
                else if(previous == term){
                    prevcoef = prevcoef + coef;
                }
                else{
                    if(prevcoef.numer!=0){
                        if(prevcoef.numer!=1){
                            if(prevcoef.numer>0)
                                out << "+";
                            if(prevcoef.numer!=-1)
                                out << prevcoef.numer;
                            else
                                out << "-";
                            if(prevcoef.denom!=1)
                                out << "1/" << prevcoef.denom;
                            else{
                                if(prevcoef.numer!=1 && prevcoef.numer!=-1){
                                    out << "*";
                                }
                            }

                        }
                        else{
                            if(prevcoef.denom!=1)
                                out << "+1/" << prevcoef.denom << "*";
                            else
                                out << "+";
                        }
                        out << previous;
                    }

                    previous = term;
                    prevcoef = coef;
                }
            }
        }
	}
	string output = out.str();
    if(output.substr(0,2)=="+0"||output.substr(0,2)=="-0")
        output.erase(0,2);
    if(output[0]=='0')
        output.erase(0,1);
    if(output[0]=='+')
        output.erase(0,1);
    return output;
}

string Canonizer::multiplyByMinusOne(string s)
{
  //  cout << "Called multiplyByMinusOne on input '" << s << "'" << endl;
    for(int j=s.length()-1;j>=0;j--){
        if(isComparison(s[j])){
            if(s[j+1]!='+' && s[j+1]!='-'){
                s = s.substr(0,j+1) + "+" + s.substr(j+1);
            }
            break;
        }
    }
    int par = 0, len=s.length();
    for(int i=0;i<len;i++){
        if(s[i]=='(')
            par++;
        if(s[i]==')')
            par--;
        if(par==0 && s[i]=='-')
            s[i] = '+';
        else if(par==0 && s[i]=='+')
            s[i] = '-';
    }
    if(s[0]=='+')
        s.erase(0,1);
    while(!isComparison(s[len]))
        len--;
    if(s[len]=='>')
        s[len] = '<';
    else if(s[len]=='<')
        s[len] = '>';
    if(s[len-1]=='>')
        s[len-1] = '<';
    else if(s[len-1]=='<')
        s[len-1] = '>';
    if(s[len+1]=='+')
        s.erase(len+1,1);

    if(s[s.length()-2]=='-' && s.back()=='0'){
        s.erase(s.length()-2,1);
    }


   // cout << "returning: " << s << endl;
    return s;
}



string Canonizer::takeToRHS(string s)
{
  //  cout << "Called takeToRHS on input '" << s << "'" << endl;
    if(s[0]=='-')
        s = multiplyByMinusOne(s);
    if(s[0]=='+')
        s.erase(0,1);
    //if(islower(s[0]))
    //    return s;
    int i = 0;
    while(s[i]!='+' && s[i]!='-' && !isComparison(s[i])){
        i++;
    }

   // cout << "invids.at(" << takevarnum(s,0) << ") = " << invIDs.at(takevarnum(s,0)) << endl;

    if( (s[0]=='X'||s[0]=='V') && !isVariable(invIDs.at(takevarnum(s,0))) ){
       // cout << "takeToRHS stopping as invids.at(" << takevarnum(s,0) << ") = " << invIDs.at(takevarnum(s,0)) << " is not an ASP var" << endl;

        return s;
    }

    stringstream rhs;
    rhs << "-" << s.substr(0,i);
    s.erase(0,i);



    if(isComparison(s[s.length()-2]) && s.back()=='0')
        s.pop_back();

    s.append(rhs.str());

    if(isComparison(s[0])){
        s = "0" + s;
       // cout << "returning RHS-only '" << s << "'" << endl;
        return s;
    }


    return takeToRHS(s);
}

string Canonizer::canonicalForm(string s, set<int> &scope)
{


    s = makeLHSonly(s);
 //   cout << "\nmakeLHSonly -> " << s << endl;

    string rhs = "";
    for(int i=s.length()-2;i>=0;i--){
        if(!isComparison(s[i])){
            rhs = s.substr(i+1);
            s.erase(i+1);
            break;
        }
    }
    //cout << "RHS stored apart: " << rhs << endl;

    s = makeVariablesIDs(s,scope);
   // cout << "\nmakeVariablesIDs -> " << s << endl;
    s = sortExpr(s);
 //   cout << "\nsortExpr -> " << s << endl;
    s = computeSums(s);
  //  cout << "\ncomputeSums -> " << s << endl;

    s.append(rhs);
  //  cout << "\nAppend zero RHS -> " << s << endl;
    s = takeToRHS(s);
  //  cout << "\ntakeToRHS -> " << s << endl;

    return s;
}

string Canonizer::canonicalForm(string s)
{
  //  cout << "Called canonicalForm (dummy scope) on input '" << s << "'" << endl;

    s = makeLHSonly(s);
 //   cout << "\nmakeLHSonly -> " << s << endl;

    string rhs = "";
    for(int i=s.length()-2;i>=0;i--){
        if(!isComparison(s[i])){
            rhs = s.substr(i+1);
            s.erase(i+1);
            break;
        }
    }
  //  cout << "RHS stored apart: " << rhs << endl;
    set<int> dummy_scope;
    s = makeVariablesIDs(s,dummy_scope);
  //  cout << "\nmakeVariablesIDs -> " << s << endl;
    s = sortExpr(s);
  //  cout << "\nsortExpr -> " << s << endl;
    s = computeSums(s);
  //  cout << "\ncomputeSums -> " << s << endl;

    s.append(rhs);
  //  cout << "\nAppend zero RHS -> " << s << endl;
    s = takeToRHS(s);
 //   cout << "\ntakeToRHS -> " << s << endl;
    s = backToVarNames(s);
  //  cout << "\nbackToVarNames -> " << s << endl;
    boost::replace_all(s,"(","{");
    boost::replace_all(s,")","}");
    boost::replace_all(s,",",";");
    //cout << "\nRestore symbols -> " << s << endl;

    return s;
}

string Canonizer::backToVarNames(string s)
{

    stringstream out;
    bool v = false;
    for(int i=0;i<s.length();i++){
        if(s[i]=='X'){
            out << invIDs.at(takevarnum(s,i));
            v = true;
        }
        else if(s[i]=='V'){
            out << "(1/" << invIDs.at(takevarnum(s,i)) << ")";
            v = true;
        }
        else if(v && !isdigit(s[i])){
            v = false;
            out << s[i];
        }
        else if(!v){
            out << s[i];
        }
    }


    return out.str();
}

expression Canonizer::makeExpression(string in, string id)
{
    expression ex;
    ex.ID = id;
    ex.raw = in;
    ex.scope.clear();
    //cout << "Called 'makeExpression' on input in = " << in << ", id = " << id << endl;
    string canon = canonicalForm(in,ex.scope);
    for(int i=0;i<canon.length();i++){
        if(isComparison(canon[i])){
            ex.lhs = canon.substr(0,i);
            if(canon[i]=='!'){
                ex.co = NEQ;
                ex.RHS = canon.substr(i+2);
                if(isNumber(ex.RHS)){
                    ex.intrhs = true;
                    ex.rhs = atoi(ex.RHS.c_str());
                }
                break;
            }
            else if(canon[i]=='='){
                ex.co = EQ;
                ex.RHS = canon.substr(i+1);
                if(isNumber(ex.RHS)){
                    ex.intrhs = true;
                    ex.rhs = atoi(ex.RHS.c_str());
                }
                break;
            }
            else if(canon[i]=='<'){
                if(canon[i+1]=='='){
                    ex.co = LEQ;
                    ex.RHS = canon.substr(i+2);
                    if(isNumber(ex.RHS)){
                        ex.intrhs = true;
                        ex.rhs = atoi(ex.RHS.c_str());
                    }
                    break;
                }
                else{
                    ex.co = LT;
                    ex.RHS = canon.substr(i+1);
                    if(isNumber(ex.RHS)){
                        ex.intrhs = true;
                        ex.rhs = atoi(ex.RHS.c_str());
                    }
                    break;
                }
            }
            else{
                if(canon[i+1]=='='){
                    ex.co = GEQ;
                    ex.RHS = canon.substr(i+2);
                    if(isNumber(ex.RHS)){
                        ex.intrhs = true;
                        ex.rhs = atoi(ex.RHS.c_str());
                    }
                    break;
                }
                else{
                    ex.co = GT;
                    ex.RHS = canon.substr(i+1);
                    if(isNumber(ex.RHS)){
                        ex.intrhs = true;
                        ex.rhs = atoi(ex.RHS.c_str());
                    }
                    break;
                }
            }
        }
    }

    return ex;
}

string Canonizer::generalizeNongroundExpr(string s)
{
    char_separator<char> sep(" ", "+-*/{};=<>!", drop_empty_tokens);
	tokenizer<char_separator<char> > tokens(s, sep);

	stringstream out;
	int n = 0;
	for ( tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
        string token = *it;
        if(isVariable(token)){
            out << "Var" << n;
            n++;
        }
        else{
            out << token;
        }
	}

	return out.str();
}

string Canonizer::generalizeVars(string s, int arity)
{
    for(int i=0;i<arity;i++){
        stringstream v;
        v << "Var" << i;
        boost::replace_all(s,v.str(),"?");
        v.str(string());
    }
    return s;
}

void Canonizer::addExpression(string s, string id)
{
    Expressions.push_back(makeExpression(s,id));
}

string Canonizer::analize()
{
    //cout << "Called 'Canonizer::ground_analyze'" << endl;
    if(!Expressions.back().intrhs)
        return "\nERROR: nonground expression found by ground_analize\n";

    for(int i=0;i<Expressions.size()-1;i++){
        if(!Expressions[i].intrhs)
            return "\nERROR: nonground expression found by ground_analize\n";
        if(Expressions.back().lhs == Expressions[i].lhs){
            if(Expressions.back().co == Expressions[i].co){
                if(Expressions.back().rhs == Expressions[i].rhs){
                    //cout << Expressions.back().ID << " EQUIVALENT " << Expressions[i].ID << endl;

                    stringstream out;

                    out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << ".\n"
                        << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;

                    Expressions.pop_back();
                    return out.str();
                }
                else if(Expressions.back().rhs < Expressions[i].rhs){
                    if(Expressions.back().co == EQ){
                        //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                        stringstream out;

                        out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;

                        return out.str();

                    }
                    if(Expressions.back().co == LT || Expressions.back().co == LEQ){
                        //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;
                        stringstream out;
                        out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                        Expressions[i] = Expressions.back();    Expressions.pop_back();
                        return out.str();
                    }
                    if(Expressions.back().co == GT || Expressions.back().co == GEQ){
                        //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;
                        stringstream out;
                        out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                        Expressions.pop_back();
                        return out.str();
                    }
                }
                else if(Expressions.back().rhs > Expressions[i].rhs){
                    if(Expressions.back().co == EQ){
                        //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                        stringstream out;
                        out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                        return out.str();
                    }
                    if(Expressions.back().co == LT || Expressions.back().co == LEQ){
                        //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                        stringstream out;
                        out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                        Expressions.pop_back();
                        return out.str();
                    }
                    if(Expressions.back().co == GT || Expressions.back().co == GEQ){
                        //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                        stringstream out;
                        out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                        Expressions[i] = Expressions.back();    Expressions.pop_back();
                        return out.str();
                    }
                }
            }
            else{
                if(Expressions.back().rhs == Expressions[i].rhs){
                    if(Expressions.back().co == EQ){
                        if(Expressions[i].co == GEQ || Expressions[i].co == LEQ){
                             //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();
                        }
                    }
                    if(Expressions.back().co == NEQ){
                        if(Expressions[i].co == GT || Expressions[i].co == LT){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();

                        }
                    }
                    if(Expressions.back().co == LT){
                        if(Expressions[i].co == LEQ || Expressions[i].co == NEQ){
                             //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();
                        }
                    }
                    if(Expressions.back().co == GT){
                        if(Expressions[i].co == GEQ || Expressions[i].co == NEQ){
                             //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();

                        }
                    }
                    if(Expressions.back().co == LEQ){
                        if(Expressions[i].co == LT || Expressions[i].co == EQ){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();

                        }
                    }
                    if(Expressions.back().co == GEQ){
                        if(Expressions[i].co == GT || Expressions[i].co == EQ){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();

                        }
                    }
                }
                else if(Expressions.back().rhs < Expressions[i].rhs){
                    if(Expressions.back().co == EQ || Expressions.back().co == LEQ || Expressions.back().co == LT){
                        if(Expressions[i].co == LT || Expressions[i].co == LEQ || Expressions[i].co == NEQ){
                            //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();

                        }
                    }
                    if(Expressions.back().co == GEQ || Expressions.back().co == GT || Expressions.back().co == NEQ){
                        if(Expressions[i].co == GT || Expressions[i].co == GEQ || Expressions[i].co == EQ){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " COMPATIBLE " << Expressions[i].ID;

                            relate(Expressions.back(),Expressions[i]);
                            relate(Expressions[i],Expressions.back());

                            if(Expressions.back().hasmin && Expressions.back().hasmax && chances(Expressions.back())==0){
                                stringstream out;
                                out << ":- ";
                                for(int j=0;j<Expressions.size()-1;j++)
                                    if(Expressions[j].lhs == Expressions.back().lhs)
                                        out << Expressions[j].ID << ", ";
                                out << Expressions.back().ID << "." << endl;
                                return out.str();
                            }
                        }
                    }
                }
                else{
                    if(Expressions.back().co == EQ || Expressions.back().co == GEQ || Expressions.back().co == GT){
                        if(Expressions[i].co == GT || Expressions[i].co == GEQ || Expressions[i].co == NEQ){
                            //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();

                        }
                    }
                    if(Expressions.back().co == LEQ || Expressions.back().co == LT || Expressions.back().co == NEQ){
                        if(Expressions[i].co == LT || Expressions[i].co == LEQ || Expressions[i].co == EQ){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                        else{
                            //cout << Expressions.back().ID << " COMPATIBLE " << Expressions[i].ID;

                            relate(Expressions.back(),Expressions[i]);
                            relate(Expressions[i],Expressions.back());

                            if(Expressions.back().hasmin && Expressions.back().hasmax && chances(Expressions.back())==0){
                                stringstream out;
                                out << ":- ";
                                for(int j=0;j<Expressions.size()-1;j++)
                                    if(Expressions[j].lhs == Expressions.back().lhs)
                                        out << Expressions[j].ID << ",";
                                out << Expressions.back().ID << "." << endl;
                                return out.str();
                            }
                        }
                    }
                }
            }
        }
        /*else{
            cout << Expressions.back().ID << " DIFFERENT_LHS " << Expressions[i].ID << endl;
        }*/
    }

    return "";

//    for(int i=0;i<Expressions.size();i++){
//        cout << "Expressions[" << i << "] :" << endl
//             << "ID = " << Expressions[i].ID << endl
//             << "raw = " << Expressions[i].raw << endl
//             << "lhs = " << Expressions[i].lhs << endl
//             << "co = ";
//             if(Expressions[i].co==EQ)
//                cout << "=" << endl;
//             else if(Expressions[i].co==NEQ)
//                cout << "!=" << endl;
//             else if(Expressions[i].co==LT)
//                cout << "<" << endl;
//             else if(Expressions[i].co==LEQ)
//                cout << "<=" << endl;
//             else if(Expressions[i].co==GEQ)
//                cout << ">=" << endl;
//             else
//                cout << ">" << endl;
//
//             cout << "rhs = " << Expressions[i].rhs << endl
//             << "scope = {";
//        for(auto it=Expressions[i].scope.begin();it!=Expressions[i].scope.end();it++)
//                cout << *it << " ";
//        cout << "}" << endl << "\n" << endl;
//    }


}

int Canonizer::computeMinVal(expression a)
{
    string e = a.lhs;
    bool Plus = true;
    double res = 0;
    boost::char_separator(" ","+-",boost::drop_empty_tokens) sep;
    boost::tokenizer<boost::char_separator<char> > tokens(e, sep);
    for ( boost::tokenizer<boost::char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
        string tok = *it;
        if(tok=="+")
            Plus = true;
        if(tok=="-")
            Plus = false;
        else{
            bool Times = true;
            double coef = 1;
            vector<int> scp;
            boost::char_separator(" ","*/",boost::drop_empty_tokens) sep1;
            boost::tokenizer<boost::char_separator<char> > tokens1(tok, sep1);
            for ( boost::tokenizer<boost::char_separator<char> >::iterator jt = tokens1.begin(); jt != tokens1.end(); ++jt) {
                string tok1 = *jt;
                if(tok1=="*")
                    Times = true;
                if(tok1=="/")
                    Times = false;
                else{
                    if(isNumber(tok1)){
                        if(Times)
                            coef *= atof(tok1.c_str());
                        else{
                            coef /= atof(tok1.c_str());
                            Times = true;
                        }
                    }
                    else{
                        if(tok1[1]=='X'){
                            if(Plus){
                                res += coef * dom_min;
                            }
                            else{
                                res += coef * dom_max;
                            }
                        }
                        else if(tok1[0]=='V'){
                            if(Plus){
                                res += coef * ( 1 / dom_max );
                            }
                            else{
                                res += coef * ( 1 / dom_min);
                            }
                        }

                    }
                }
            }
        }
    }
    return floor(res)+1;

}

/*string Canonizer::renameRHSvars(string RHS, int &n)
{
    map<int,int> varmap;
    for(int i=2;i<RHS.length()-1;i++){
        if(RHS.substr(i-2,3)=="Var"){
            int m = takevarnum(RHS,i);
            if(varmap.count(m)==0){
                varmap.insert({m,n});
                n++;
            }
            stringstream t;
            t << "N" << varmap.at(m);
            RHS = RHS.substr(0,i-2) + t.str() + RHS.substr(i+numofdgts(varmap.at(m))+1);
        }
    }
    return RHS;
}*/


/*expression Canonizer::makeExprForNewRule(expression a, int &nv)
{
    a.RHS = renameRHSvars(a.RHS,nv);

    int i = 8;
    while(!isComparison(a.ID[i]))
        i++;
    while(isComparison(a.ID[i]))
        i++;
    size_t qe = a.ID.rfind('"');
    a.ID = a.ID.substr(0,i) + a.RHS + "\"";
    stringstream t;
    for(int j=0;j<nv;j++){
        t << ",N" << j;
        a.ID.append(t.str());
        t.str(string());
    }
    a.ID.append(")");

    return a;
}*/

/*string Canonizer::nonground_analize()
{
    //cout << "Called 'Canonizer::nonground_analyze'" << endl;

    for(int i=0;i<Expressions.size()-1;i++){
        if(Expressions.back().intrhs && Expressions[i].intrhs && ){
            if(Expressions.back().lhs == Expressions[i].lhs){
                if(Expressions.back().co == Expressions[i].co){
                    if(Expressions.back().rhs == Expressions[i].rhs){
                        //cout << Expressions.back().ID << " EQUIVALENT " << Expressions[i].ID << endl;

                        stringstream out;

                        out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << ".\n"
                            << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;

                        Expressions.pop_back();
                        return out.str();
                    }
                    else if(Expressions.back().rhs < Expressions[i].rhs){
                        if(Expressions.back().co == EQ){
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;

                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;

                            return out.str();

                        }
                        if(Expressions.back().co == LT || Expressions.back().co == LEQ){
                            //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;
                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                        if(Expressions.back().co == GT || Expressions.back().co == GEQ){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;
                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                    }
                    else if(Expressions.back().rhs > Expressions[i].rhs){
                        if(Expressions.back().co == EQ){
                            //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                            return out.str();
                        }
                        if(Expressions.back().co == LT || Expressions.back().co == LEQ){
                            //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                            Expressions.pop_back();
                            return out.str();
                        }
                        if(Expressions.back().co == GT || Expressions.back().co == GEQ){
                            //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                            stringstream out;
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                            Expressions[i] = Expressions.back();    Expressions.pop_back();
                            return out.str();
                        }
                    }
                }
                else{
                    if(Expressions.back().rhs == Expressions[i].rhs){
                        if(Expressions.back().co == EQ){
                            if(Expressions[i].co == GEQ || Expressions[i].co == LEQ){
                                 //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                                Expressions[i] = Expressions.back();    Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();
                            }
                        }
                        if(Expressions.back().co == NEQ){
                            if(Expressions[i].co == GT || Expressions[i].co == LT){
                                //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                                Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();

                            }
                        }
                        if(Expressions.back().co == LT){
                            if(Expressions[i].co == LEQ || Expressions[i].co == NEQ){
                                 //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                                Expressions[i] = Expressions.back();    Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();
                            }
                        }
                        if(Expressions.back().co == GT){
                            if(Expressions[i].co == GEQ || Expressions[i].co == NEQ){
                                 //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                                Expressions[i] = Expressions.back();    Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();

                            }
                        }
                        if(Expressions.back().co == LEQ){
                            if(Expressions[i].co == LT || Expressions[i].co == EQ){
                                //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                                Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();

                            }
                        }
                        if(Expressions.back().co == GEQ){
                            if(Expressions[i].co == GT || Expressions[i].co == EQ){
                                //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                                Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();

                            }
                        }
                    }
                    else if(Expressions.back().rhs < Expressions[i].rhs){
                        if(Expressions.back().co == EQ || Expressions.back().co == LEQ || Expressions.back().co == LT){
                            if(Expressions[i].co == LT || Expressions[i].co == LEQ || Expressions[i].co == NEQ){
                                //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                                Expressions[i] = Expressions.back();    Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();

                            }
                        }
                        if(Expressions.back().co == GEQ || Expressions.back().co == GT || Expressions.back().co == NEQ){
                            if(Expressions[i].co == GT || Expressions[i].co == GEQ || Expressions[i].co == EQ){
                                //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                                Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " COMPATIBLE " << Expressions[i].ID;

                                relate(Expressions.back(),Expressions[i]);
                                relate(Expressions[i],Expressions.back());

                                if(Expressions.back().hasmin && Expressions.back().hasmax && chances(Expressions.back())==0){
                                    stringstream out;
                                    out << ":- ";
                                    for(int j=0;j<Expressions.size()-1;j++)
                                        if(Expressions[j].lhs == Expressions.back().lhs)
                                            out << Expressions[j].ID << ", ";
                                    out << Expressions.back().ID << "." << endl;
                                    return out.str();
                                }
                            }
                        }
                    }
                    else{
                        if(Expressions.back().co == EQ || Expressions.back().co == GEQ || Expressions.back().co == GT){
                            if(Expressions[i].co == GT || Expressions[i].co == GEQ || Expressions[i].co == NEQ){
                                //cout << Expressions.back().ID << " IMPLIES " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << "." << endl;
                                Expressions[i] = Expressions.back();    Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " INCOMPATIBLE " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions.back().ID << ", " << Expressions[i].ID << "." << endl;
                                return out.str();

                            }
                        }
                        if(Expressions.back().co == LEQ || Expressions.back().co == LT || Expressions.back().co == NEQ){
                            if(Expressions[i].co == LT || Expressions[i].co == LEQ || Expressions[i].co == EQ){
                                //cout << Expressions.back().ID << " IMPLIED_BY " << Expressions[i].ID << endl;

                                stringstream out;
                                out << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;
                                Expressions.pop_back();
                                return out.str();
                            }
                            else{
                                //cout << Expressions.back().ID << " COMPATIBLE " << Expressions[i].ID;

                                relate(Expressions.back(),Expressions[i]);
                                relate(Expressions[i],Expressions.back());

                                if(Expressions.back().hasmin && Expressions.back().hasmax && chances(Expressions.back())==0){
                                    stringstream out;
                                    out << ":- ";
                                    for(int j=0;j<Expressions.size()-1;j++)
                                        if(Expressions[j].lhs == Expressions.back().lhs)
                                            out << Expressions[j].ID << ",";
                                    out << Expressions.back().ID << "." << endl;
                                    return out.str();
                                }
                            }
                        }
                    }
                }
            }
            /*else{
            cout << Expressions.back().ID << " DIFFERENT_LHS " << Expressions[i].ID << endl;
            }*/
    /*    }
        else{
            if(Expressions.back().lhs == Expressions[i].lhs){
                if(Expressions.back().co == Expressions[i].co){
                    if(Expressions.back().RHS == Expressions[i].RHS){
                        stringstream out;
                        if(Expressions.back().ID != Expressions[i].ID)
                            out << ":- " << Expressions.back().ID << ", not " << Expressions[i].ID << ".\n"
                                << ":- " << Expressions[i].ID << ", not " << Expressions.back().ID << "." << endl;

                        Expressions.pop_back();
                        ///return out.str();
                    }
                    else{
                        int vvnn = 0;
                            a = makeExprForNewRule(Expressions.back(),vvnn);
                            b = makeExprForNewRule(Expressions[i],vvnn);

                        if(Expressions.back().co == EQ){
                            stringstream out;
                            out << ":- " << a.ID << "," << b.ID << "," << a.RHS << " != " << b.RHS << "." << endl;
                            ///return out.str();
                        }
                        else if(Expressions.back().co == LEQ || Expressions.back().co == LT){
                            stringstream out;
                            out << ":- " << a.ID << ", not " << b.ID << "," << a.RHS << " <= " << b.RHS << ".\n"
                                << ":- " << b.ID << ", not " << a.ID << "," << b.RHS << " <= " << a.RHS << "." << endl;
                            ///return out.str();
                        }
                        else if(Expressions.back().co == GEQ || Expressions.back().co == GT){
                            stringstream out;
                            out << ":- " << a.ID << ", not " << b.ID << "," << a.RHS << " >= " << b.RHS << ".\n"
                                << ":- " << b.ID << ", not " << a.ID << "," << b.RHS << " >= " << a.RHS << "." << endl;
                            ///return out.str();
                        }
                    }
                }
                else{

                }
            }
        }

    }

    return "";


} */

void Canonizer::printExpressions()
{
    cout << "Expressions contains: " << endl;
    for(int k=0;k<Expressions.size();k++){
        cout << Expressions[k].lhs;
        if(Expressions[k].co == EQ)
            cout << "=";
        else if(Expressions[k].co == NEQ)
            cout << "!=";
        else if(Expressions[k].co == LEQ)
            cout << "<=";
        else if(Expressions[k].co == GEQ)
            cout << ">=";
        else if(Expressions[k].co == LT)
            cout << "<";
        else
            cout << ">";
        cout << Expressions[k].rhs << endl;
    }
    cout << "----------------------\n" << endl;
}

void Canonizer::printvarIDs()
{
    cout << "VarIDs contains: " << endl;
    for(auto it=varIDs.begin();it!=varIDs.end();it++)
        cout << it->first << "<-->" << it->second << endl;
    cout << "----------------------\n" << endl;
}

void Canonizer::clear()
{
    varIDs.clear();
    invIDs.clear();
    freeASPvarnum = 0;
    Expressions.clear();
}
