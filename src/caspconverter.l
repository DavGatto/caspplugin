%option noyywrap
%option c++
%option outfile="caspconverter.cpp"


%x QUOTEDSTRING
%x COMMENT

IFSYMBOL :\-
DISJ " v "
COMPARISON (!=|<>|<=|>=|<|>|=)
BINOPERAT [\+\-\*\/]
DCOMPARISON \$(!=|<>|<=|>=|<|>|=)
DBINOPERAT \$[\+\-\*\/]
TERM [A-Za-z][0-9a-zA-Z\-\_]*
ARGUMENT \([^\.\:\|\)]+[^\.\:\|\(]*\)
ATMARG {TERM}{ARGUMENT}?
ATOM {ATMARG}(\(({ATMARG}\,)*{ATMARG}\))?
DHCOMMENT \%[.\n]*\n
ANYSPACE [ \t]*
AOD {ATOM}|[0-9]+
NODOLLAR (({BINOPERAT}{ANYSPACE})?{AOD}{ANYSPACE}{BINOPERAT}{ANYSPACE})*

%{
#include "casp/caspconverter.h"
#include "casp/utility.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <exception>

#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

using namespace std;

bool syntaxCheckCnsAtm(string s, int LineNum)
{
    int p = 0, c = 0;
//    if(s[0]=='$'||isComparison(s[0]||isBinOperator(s[0]))){
//        cerr  << "Line " << LineNum << ": In constraint atom '" << s << "', first symbol is '" << s[0] << "'. Linear expressions must not start with an operator, please fix it (maybe begin with a zero is a good idea).";
//
//        return 1;
//    }
    for(int i=0;i<s.length();i++){
        if(!isValidConstraintSym(s[i])) { cerr << "Line " << LineNum << ": Invalid symbol '" << s[i] << "' in constraint atom '" << s << "'.";  return 1;}
        if(s[i]=='{') { p++; }
        if(s[i]=='}') { p--; if(p<0){ cerr << "Line " << LineNum << ": Unbalanced parenthesis in constraint atom'" << s << "'.";  return 1; }}
        if(s[i]=='$' && ( !isComparison(s[i+1]) && !isBinOperator(s[i+1]) || i+1==s.length() )){
            cerr << "Line " << LineNum << ": Invalid operator symbol '" << s[i] << "' in constraint atom '" << s << "'.";
            return 1;
        }
        if( isComparison( s[i] ) ){
            c++;
            if(c>2){
                cerr << "Line " << LineNum << ": Constraint atom '" << s << "' has more than two comparison symbols.";

                return 1;
            }
            if(c==2 && (s[i]!='='||!isComparison(s[i-1]))){
                cerr << "Line " << LineNum << ": Invalid comparison symbol '" << s[i-1] << s[i] << "' in constraint atom '" << s << "'.";

                return 1;
            }
            if(s[i]=='!' && ( i+1==s.length() || s[i+1]!='=')){
                cerr << "Line " << LineNum << ": Invalid comparison symbol '" << s[i] << s[i+1] << "' in constraint atom '" << s << "'.";

                return 1;
            }
        }
    }
    if(isBinOperator(s.back()) || isComparison(s.back()) || s.back()=='$'){
        cerr <<"Line " << LineNum << ": Constraint atom '" << s << "' ends with symbol '" << s.back() <<"'.";

        return 1;
    }
    if(c<1){
        cerr << "Line " << LineNum << ": Missing comparison operator in constraint atom '" << s << "'.";

        return 1;
    }
    if(p!=0){
        cerr << "Line " << LineNum << ": Unbalanced parenthesis in constraint atom '" << s << "'.";

        return 1;
    }

    return 0;
}

int CnsNum = 0, LineNum = 1, VarNum = 0, iii = 0, len = 0, dcn = 0, Ref = 0, arity = 0;
bool flag = 0; /*body = 0;*/
char * tmp;
string info = "", TMPS = "";
stringstream out;
stringstream expr_id;
stringstream expr_str;

unordered_map<int,vector<string>> depmap = {};
unordered_map<string,int> varids = {};
unordered_map<string,int> cnsids = {};
unordered_map<int,int> arityof = {};

Canonizer dummycanon;


%}

%%

\% { BEGIN(COMMENT); ECHO; }

<COMMENT>\n { ++LineNum; BEGIN(INITIAL);}




\" { BEGIN(QUOTEDSTRING); ECHO; }

<QUOTEDSTRING>\" { BEGIN(INITIAL);}



<INITIAL,QUOTEDSTRING>\n { ++LineNum; ECHO; }

\$("domain"|"maximize"|"minimize")\(.+\) { /**yyout << "DOM FOUND: " << YYText() << endl;*/
                                            for(int k=1;k<YYLeng();k++){ if(YYText()[k]==')'){ *yyout << "\"";} *yyout << YYText()[k]; if(YYText()[k]=='('){ *yyout << "\"";}}}

({NODOLLAR}{AOD}{ANYSPACE}{DBINOPERAT}{ANYSPACE})*{NODOLLAR}{AOD}{ANYSPACE}{DCOMPARISON}{ANYSPACE}({NODOLLAR}{AOD}{ANYSPACE}{DBINOPERAT}{ANYSPACE})*{NODOLLAR}{AOD} {

          //cout << endl << "CONSTRAINT FOUND: |" << YYText() << "|" << endl;

        dcn = 0; arity = 0; out.str(string()); expr_id.str(string()); expr_str.str(string());
        tmp = new char [YYLeng()]; tmp = strdup(YYText()); len = 0;
        for(int k=0;k<YYLeng();k++){
            if(tmp[k] != ' ' && tmp[k] != '\t' && tmp[k] != '\0'){ len++; }
            if(tmp[k] == '$'){ dcn++;}
        }
        for(int k=0;k<YYLeng();k++){
            /** Change this section accordingly with the fashion you prefer for the constraint string in output: */
            if(tmp[k] != ' ' && tmp[k] != '\t'){
                /*if(tmp[k]=='$'){out << ' ';}*/
                /*else{out << tmp[k];}*/ /** Print space in place of $ */
                if(tmp[k] == '$'){ out << '$';}  /** Print null char in place of $ */
                else if(tmp[k] == '('){ out << '{';}
                else if(tmp[k] == ')'){ out << '}';}
                else if(tmp[k] == ','){ out << ';';}
                else { out << tmp[k];}
                len--;
                if(len==0){break;}
            }
        }
         /** Syntax check for constraint atoms: */
        TMPS=out.str();

        if(syntaxCheckCnsAtm(TMPS,LineNum)){ cerr << endl << "CASPplugin syntax error" << endl; throw logic_error("Syntax error in constraint atom."); /*throw dlvhex::PluginError("\n CASPplugin syntax error! \n");*/}


        for(int k=0;k<TMPS.length();k++){ if(TMPS[k]=='$'){ TMPS.erase(TMPS.begin()+k); } }

        dummycanon.clear();
        TMPS = dummycanon.canonicalForm(TMPS);
        info = TMPS;
        //cout << "TMPS in canonical form: " << TMPS << endl;
        //info = dummycanon.generalizeNongroundExpr(TMPS);
        //cout << "\tinfo = " << info << endl;

        if(cnsids.count(info)==0){ CnsNum++; cnsids.insert({info,CnsNum}); arityof.insert({CnsNum,0}); vector<string> vars; depmap.insert({CnsNum,vars});}
        Ref = cnsids[info];

        expr_id << "expr_" << Ref << "_";
        expr_str << "(\"" << TMPS << "\"";


        boost::char_separator<char> sep("+-*/=<>! ", "", boost::drop_empty_tokens);
        boost::tokenizer<boost::char_separator<char> > tokens(TMPS, sep);

        for ( boost::tokenizer<boost::char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
            string tok = *it;
            boost::replace_all(tok,"{","(");
            boost::replace_all(tok,"}",")");
            boost::replace_all(tok,";",",");
            if(!isNumberOrPar(tok)){
                //cout << endl << "VARIABLE FOUND: " << tok << endl;
                expr_str << "," << extractVariables(tok,arity); /*nguess << extractVariables(out.str());*/
                string tokvars = generalizeVariables(tok); /*cout << endl << TMPS << endl;*/
                if(tokvars!="?"){
                    depmap.at(Ref).push_back(tokvars);
                    if(varids.count(tokvars) == 0){
                        varids.insert({tokvars,VarNum}); VarNum++;
                    }
                }
            }
        }
        expr_str << ")";
        //expr_id << arity;

        *yyout << expr_id.str() << arity << expr_str.str();
        //cout << "CONVERTED: " << expr_id.str() << arity << expr_str.str();
        //cout << "Wrote: " << expr_id.str() << expr_str.str();
        //int foo=0;
        //string genvs = extractVariables(info,foo);
        //genvs = expr_id.str() + "(\"" + info + "\"," + genvs + ")";
        //cout << "genvs = " << genvs << endl;
        //canon.addExpression(info,genvs);

        //canon.printvarIDs();

        arityof.at(Ref) = arity;

        delete tmp; out.str(string());
    }


%%

DefaultConverter::DefaultConverter() {
}

DefaultConverter::~DefaultConverter() {
}

void DefaultConverter::convert(istream& i, ostream& o) {
	string input;
	while (getline(i, input)) {
		o << input << endl;
	}
}

CaspConverter::CaspConverter() {

}

CaspConverter::~CaspConverter() {

}



void CaspConverter::convert(istream& i, ostream& o)
{
    cout << "\nHello! Using modified caspconverter..." << endl;

    istream * inp = &i; ostream * outp = &o;


    /** Flex scanner call: */
    FlexLexer* lexer = new yyFlexLexer;
    while(lexer->yylex(inp, outp) != 0);


/*    canon.printExpressions();

    o << "% Nonground nogoods derived by converter:" << endl;
    o << canon.analize() << endl << "% End of nogoods from converter." << endl;

    canon.printExpressions();
    canon.printvarIDs();*/

    /** Build the constraint/variables dependency map 'Map': */
    unordered_map<int,vector<bool>> Map = {};
    for (int h=1; h<depmap.size()+1; h++){ vector<bool> vbool; Map.insert({h,vbool}); for (int j=0;j<VarNum;j++){ Map.at(h)[j] = 0;}}

    for(int i=1;i<=CnsNum;i++){
        for(int j=0;j<depmap[i].size();j++){
            Map.at(i)[varids.at(depmap.at(i)[j])] = 1;
        }
    }

    /** Uncomment to print Map: */
//    cout << "\n\nMatrice delle adiacenze:" << endl;
//    for ( auto it = Map.begin(); it != Map.end(); ++it ){
//        cout << "{" << it->first << ": ";
//        for (int j=0;j<VarNum;j++){
//            cout << it->second[j] << " ";
//        }
//        cout << "}" << endl;
//    }

    /** Uncomment to print varids (Variables IDentificators): */
//    cout << "\nVarIDs:" << endl;
//    for ( auto it = varids.begin(); it != varids.end(); ++it ){
//        cout << it->second << ": " << it->first << endl;
//    }

    /** Constraints' interdependence analysis: */
    bool erased = 1, matched = 0;
    vector<int> list; len = 1;
    while( Map.size() != 0 ){
        while( Map.count(len) == 0 ){ len++; }
        erased = 0;
        if( Map.size() == 1 ){
            list.push_back(len);
            Map.erase(len);
            break;
        }
        for(int i=len+1; i<CnsNum+1; i++){
            matched = 0;
            if( Map.count(i) != 0){
                for(int j=0; j<VarNum; j++){
                    if(matched) {break;}
                    if( Map.at(len)[j] && Map.at(i)[j] ){
                        matched = 1;
                        for(int k=j+1; k<VarNum; k++){
                            if(Map.at(i)[k]){
                                Map.at(len)[k] = 1;
                                Map.at(i)[k] = 0;
                            }
                        }
                    }
                }
                if(matched){
                    list.push_back(i);
                    Map.erase(i);
                    erased = 1;
                }
            }
        }
        if(!erased){
            list.push_back(len);
            list.push_back(0);
            Map.erase(len);
            len++;
        }
    }
    /** Uncomment to print Map (which now should be empty):
    cout << "\n\nMatrice delle adiacenze:" << endl;
    for ( auto it = Map.begin(); it != Map.end(); ++it ){
        cout << "{" << it->first << ": ";
        for (int j=0;j<VarNum;j++){
            cout << it->second[j] << " ";
        }
        cout << "}" << endl;
    }
    /** Uncomment to print the the vector 'list' (groups of constraints separed by zeroes, each group is independent of the other):
    cout << "\n\nLIST = ";for(int h=0; h<list.size(); h++){cout << list[h] << " ";} cout << endl;*/
    /*cout << "\n\nArities:" << endl;
    for ( auto it = arityof.begin(); it != arityof.end(); ++it ){
        cout << it->first << ": " << it->second << endl;
    }*/

    /** Print the hex calls on the output stream: */
    iii=0;
    while(iii<list.size()){
        o << "\n:- not &casp[domain,maximize,minimize";
        while(list[iii]!=0 && iii<list.size()){
            *outp << ",expr_" << list[iii] << "_" << arityof.at(list[iii]) << ",not_expr_" << list[iii] << "_" << arityof.at(list[iii]);
            iii++;
        }
        for (int h=0; h<15; h++){ o << ",dummy_sum_predicate" << h;}
        o << "]()." << endl;
        iii++;
    }
    /** End of 'convert' */
    return;
}
