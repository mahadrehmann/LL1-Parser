#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <algorithm>
#include <cctype> 

using namespace std;

const string EPSILON = "^";

// Global FIRST/FOLLOW and parsing table
unordered_map<string, unordered_set<string>> firstSet;
unordered_map<string, unordered_set<string>> followSet;
unordered_map<string, unordered_map<string, vector<string>>> parsingTable;


void printFirstSets();

// A rule: LHS and a list of productions.
struct Rule {
    string lhs;
    vector<string> productions;
};

// Trim whitespace from both ends of a string.
string trim(const string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if(start == string::npos)
        return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Reading the CFG from a file.
void readGrammar(const string &filename, vector<Rule> &grammar) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "File nahi Khuli" << filename << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while(getline(fin, line)) {
        line = trim(line);
        if(line.empty())
            continue;

        size_t pos = line.find("->");
        if(pos == string::npos)
            continue;

        string lhs = trim(line.substr(0, pos));
        string rhsLine = trim(line.substr(pos + 2));

        // Split the RHS by '|'
        vector<string> prods;
        size_t start = 0;
        while(true) {
            size_t pipePos = rhsLine.find('|', start);
            string prod;
            if(pipePos == string::npos) {
                prod = trim(rhsLine.substr(start));
                if(!prod.empty()) prods.push_back(prod);
                break;
            } else {
                prod = trim(rhsLine.substr(start, pipePos - start));
                if(!prod.empty()) prods.push_back(prod);
                start = pipePos + 1;
            }
        }

        // cout<<"about to push: "<<lhs <<" "<<prods[0]<<endl;

        grammar.push_back({lhs, prods});
    }
    fin.close();
}

// Print the current grammar.
void printGrammar(const vector<Rule> &grammar, ostream &print) {
    for (const auto &rule : grammar) {
        print << rule.lhs << " -> ";
        for (size_t i = 0; i < rule.productions.size(); i++) {
            print << rule.productions[i];
            if(i != rule.productions.size() - 1)
                print << " | ";
        }
        print << endl;
    }
    print << "-----------------------------------------------------\n";
}

// Removing left factoring
void leftFactoring(vector<Rule> &grammar) {
    cout << "\nPerforming Left Factoring..." << endl;
    vector<Rule> newGrammar;
    unordered_map<string, int> nonTerminalCount;

    for (auto &rule : grammar) {
        vector<string> productions = rule.productions;
        string lhs = rule.lhs;
        
        while (true) {
            string commonPrefix = "";
            bool found = false;

            for (size_t i = 0; i < productions.size(); i++) {
                for (size_t j = i + 1; j < productions.size(); j++) {
                    string prefix = "";
                    size_t minLen = min(productions[i].size(), productions[j].size());
                    for (size_t k = 0; k < minLen; k++) {
                        if (productions[i][k] == productions[j][k]) {
                            prefix += productions[i][k];
                        } else {
                            break;
                        }
                    }
                    if (!prefix.empty() && prefix.length() > commonPrefix.length()) {
                        commonPrefix = prefix;
                        found = true;
                    }
                }
            }
            
            if (!found) {
                newGrammar.push_back({lhs, productions});
                break;
            }
            
            nonTerminalCount[lhs]++;
            string newNonTerminal = lhs + to_string(nonTerminalCount[lhs]);
            vector<string> newRuleProductions;
            vector<string> remainingProductions;
            
            for (auto &prod : productions) {
                if (prod.substr(0, commonPrefix.length()) == commonPrefix) {
                    string suffix = prod.substr(commonPrefix.length());
                    if (suffix.empty()) {
                        suffix = EPSILON;
                    }
                    newRuleProductions.push_back(suffix);
                } else {
                    remainingProductions.push_back(prod);
                }
            }
            
            remainingProductions.push_back(commonPrefix + newNonTerminal);
            productions = remainingProductions;
            newGrammar.push_back({newNonTerminal, newRuleProductions});
        }
    }
    
    grammar = newGrammar;
}

// Removing both indirect and direct left recursion.
void removeLeftRecursion(vector<Rule> &grammar) {

    int numRules = grammar.size();

    // Substituting productions for checking indirect recursion.
    for (int i = 0; i < numRules; i++) {
        for (int j = 0; j < i; j++) {

            string previousNT = grammar[j].lhs;
            vector<string> newProduction;  

            // For indirect, replace value here
            for (auto &prod : grammar[i].productions) {
                
                // Replazce if same as previos
                if (prod.substr(0, previousNT.size()) == previousNT) {
                    string restProduction = prod.substr(previousNT.size());

                    for (auto &alpha : grammar[j].productions) {
                        newProduction.push_back(alpha + restProduction);
                    }
                } else {
                    newProduction.push_back(prod);
                }
            }
            grammar[i].productions = newProduction;
        }

        // Split productions into recursive and non-recursive parts.
        vector<string> recursivePart;       // recursive parts yaani A→Aα     
        vector<string> nonRecursivePart;   // non-recursive parts β’s

        for (auto &prod : grammar[i].productions) {

            if (prod.substr(0, grammar[i].lhs.size()) == grammar[i].lhs) {

                string rest = prod.substr(grammar[i].lhs.size());
                if (rest.empty())
                    rest = EPSILON;
                recursivePart.push_back(rest);
            } else {
                nonRecursivePart.push_back(prod);
            }
        }

        // If we found any recursion, transform the rule.
        if (!recursivePart.empty()) {

            string newNonTerminal = grammar[i].lhs + "'"; //new non-terminal A'
            vector<string> updatedProduction;

            for (auto &betas : nonRecursivePart) {
                if (betas == EPSILON)
                    updatedProduction.push_back(newNonTerminal); 
                else
                    updatedProduction.push_back(betas + newNonTerminal); //A → β A'
            }

            grammar[i].productions = updatedProduction;

            Rule newRule;  //Full new rule
            newRule.lhs = newNonTerminal;

            for (auto &alphas : recursivePart) {
                if (alphas == EPSILON)
                    newRule.productions.push_back(newNonTerminal);
                else
                    newRule.productions.push_back(alphas + newNonTerminal); // A' → α A' | ε
            }

            newRule.productions.push_back(EPSILON);

            grammar.push_back(newRule);

            cout << "Left recursion removed for rule " << grammar[i].lhs << endl;
        }
    }
}

// Computing first sets for all non terminals
void computeFirstSets(vector<Rule> &grammar) {
    firstSet.clear();

    for (const auto &rule : grammar) {
        firstSet[rule.lhs] = {};
    }

    bool changed;
    do {
        changed = false;

        for (const auto &rule : grammar) {
            string nt = rule.lhs;

            for (const string &prod : rule.productions) {
                size_t oldSize = firstSet[nt].size();
                vector<string> symbols;
                string recordNonTerminals;
                size_t pos = 0;
                bool check = false;

                while (pos < prod.size()) {
                    string token;
                    if (islower(prod[pos])) {
                        while (pos < prod.size() && islower(prod[pos])) {
                            token += prod[pos++];
                        }

                        symbols.push_back(token);
                        break;
                    }

                    int i = pos;

                    while (i < prod.size() && isupper(prod[i])) {
                        recordNonTerminals += prod[i++];
                    }

                    //cout << "Recorded non terminals: " << recordNonTerminals << endl;

                    if (isupper(prod[pos])) {
                        token = prod.substr(pos, 1);
                        pos++;
                        symbols.push_back(token);
                    } else if (islower(prod[pos])) {
                        while (pos < prod.size() && islower(prod[pos])) {
                            token += prod[pos++];
                        }
                        symbols.push_back(token);
                        break;
                    } else {
                        symbols.push_back(string(1, prod[pos++]));
                    }
                }

                bool allEpsilon = true;

                for (const string &symbol : symbols) {
                    // This is for lowercase letters
                    if (!isupper(symbol[0])) {
                        firstSet[nt].insert(symbol);
                        allEpsilon = false;
                        break;
                    }
                    // Non terminals 
                    else {
                        if (firstSet.find(symbol) != firstSet.end()) {
                            for (const string &ch : firstSet[symbol]) {
                                if (ch != EPSILON)
                                    firstSet[nt].insert(ch);
                            }
                            if (firstSet[symbol].find(EPSILON) == firstSet[symbol].end()) {
                                allEpsilon = false;
                                break;
                            }
                            // else {
                            //     for (int i = 1; i < recordNonTerminals.size(); ++i) {
                            //         if (recordNonTerminals[i]) {
                            //             if (firstSet[recordNonTerminals[i]].find(EPSILON) == firstSet[recordNonTerminals[i]].end()) {
                            //                 firstSet[nt].insert(ch);
                            //             }
                            //             else {
                                            
                            //             }
                            //         }
                            //         else {
                            //             break;
                            //         }
                            //     }
                            // }
                        }
                        else {
                            firstSet[nt].insert(symbol);
                        }
                    }
                }

                if (allEpsilon) {
                    firstSet[nt].insert(EPSILON);
                }

                changed |= (oldSize != firstSet[nt].size());
            }
        }
    } while (changed);
}

// Print First Sets
void printFirstSets(ostream &print) {
    for (const auto &entry : firstSet) {
        print << "First(" << entry.first << ") = { ";
        bool first = true;
        for (const string &sym : entry.second) {
            if (!first) print << ", ";
            print << sym;
            first = false;
        }
        print << " }" << endl;
    }
    print << "-----------------------------------------------------\n";
}

// Computing follow sets for all non terminals
void computeFollowSets(vector<Rule> &grammar, const string &startSymbol) {

    for (const auto &rule : grammar) {
        followSet[rule.lhs] = {};
    }
    if (!grammar.empty())
        followSet[startSymbol].insert("$");

    bool changed;
    do {
        changed = false;

        //for form A → α B β
        for (const auto &rule : grammar) {
            string A = rule.lhs;
          
            for (const string &prod : rule.productions) {

                vector<string> symbols; // For tokenizing like ["T","E'"]
                size_t pos = 0;
                while (pos < prod.size()) {
                    string token;
                    
                    //Non-terminal: uppercase
                    if (isupper(prod[pos])) {
                        token.push_back(prod[pos++]);
                        while (pos < prod.size() && (isdigit(prod[pos]) || prod[pos] == '\'')) {
                            token.push_back(prod[pos++]);
                        }
                        symbols.push_back(token);
                    }
                    // Terminal: lowercase
                    else if (islower(prod[pos])) {
                        while (pos < prod.size() && islower(prod[pos])) {
                            token.push_back(prod[pos++]);
                        }
                        symbols.push_back(token);
                    }
                    else {
                        symbols.push_back(string(1, prod[pos++]));
                    }
                }
                
                // Any non-terminal, upadate its followset
                for (size_t i = 0; i < symbols.size(); i++) {
                    if (isupper(symbols[i][0])) {
                        size_t beforeSize = followSet[symbols[i]].size();
                        
                        // first nikalne ki logic i.e  FIRST(β)
                        unordered_set<string> firstBeta;
                        bool allEpsilon = true;
                        for (size_t j = i + 1; j < symbols.size(); j++) {
                            string token = symbols[j];
                            if (!isupper(token[0])) {
                                firstBeta.insert(token);
                                allEpsilon = false;
                                break;
                            }
                            else {
                                for (const string &sym : firstSet[token]) {
                                    if (sym != EPSILON)
                                        firstBeta.insert(sym);
                                }
                                if (firstSet[token].find(EPSILON) == firstSet[token].end()) { //.end() means end tak check karo
                                    allEpsilon = false;
                                    break;
                                }
                            }
                        }
                        
                        //If no epsilon, best scenes
                        for (const string &sym : firstBeta) {
                            if (sym != EPSILON)
                                followSet[symbols[i]].insert(sym);
                        }
                        
                        //If no beta, or it has epsilon, then we add follow of A to it
                        if ((i == symbols.size() - 1) || allEpsilon) {
                            for (const string &sym : followSet[A]) {
                                followSet[symbols[i]].insert(sym);
                            }
                        }
                        if (followSet[symbols[i]].size() != beforeSize)
                            changed = true;
                    }
                }
            }
        }
    } while (changed);
}

// Print FOLLOW sets.
void printFollowSets(ostream &print) {
    for (const auto &entry : followSet) {
        print << "Follow(" << entry.first << ") = { ";
        bool first = true;
        for (const string &sym : entry.second) {
            if (!first) print << ", ";
            print << sym;
            first = false;
        }
        print << " }" << endl;
    }
    print << "-----------------------------------------------------\n";
}

// Function to construct LL(1) parsing table
unordered_map<string, unordered_map<string, vector<string>>> constructLL1ParsingTable(
    const vector<Rule>& grammar, 
    const unordered_map<string, unordered_set<string>>& firstSet, 
    const unordered_map<string, unordered_set<string>>& followSet) {
    
    unordered_map<string, unordered_map<string, vector<string>>> parsingTable;

    for (const auto &rule : grammar) {
        string nt = rule.lhs;

        for (const string &prod : rule.productions) {
            unordered_set<string> firstOfProd;
            bool containsEpsilon = false;

            size_t pos = 0;
            while (pos < prod.size()) {
                string symbol;
                
                if (islower(prod[pos])) {
                    while (pos < prod.size() && islower(prod[pos])) {
                        symbol += prod[pos++];
                    }
                } else {
                    symbol = string(1, prod[pos++]);
                }

                if (!isupper(symbol[0])) {  
                    firstOfProd.insert(symbol);
                    break;
                } else {
                    if (firstSet.count(symbol)) {
                        for (const string &ch : firstSet.at(symbol)) {
                            if (ch != EPSILON)
                                firstOfProd.insert(ch);
                        }
                        if (firstSet.at(symbol).find(EPSILON) == firstSet.at(symbol).end()) {
                            break;
                        } else {
                            containsEpsilon = true;
                        }
                    }
                }
            }

            if (containsEpsilon) {
                firstOfProd.insert(EPSILON);
            }

            for (const string& terminal : firstOfProd) {
                if (terminal != EPSILON) {
                    parsingTable[nt][terminal] = {prod};
                }
            }

            if (firstOfProd.find(EPSILON) != firstOfProd.end()) {
                const unordered_set<string>& followOfNT = followSet.at(nt);
                for (const string& terminal : followOfNT) {
                    parsingTable[nt][terminal] = {EPSILON};
                }
            }
        }
    }
    
    return parsingTable;
}
 
// Function to print the parsing table
void printParsingTable(const unordered_map<string, unordered_map<string, vector<string>>>& parsingTable, ostream &print) {
    unordered_set<string> terminals;

    for (const auto& row : parsingTable) {
        for (const auto& entry : row.second) {
            terminals.insert(entry.first);
        }
    }

    vector<string> terminalList(terminals.begin(), terminals.end());
    sort(terminalList.begin(), terminalList.end());

    print << setw(15) << left << "Non-Terminal";
    for (const string& terminal : terminalList) {
        print << setw(15) << left << terminal;
    }
    print << "\n" << string(15 + 15 * terminalList.size(), '-') << "\n";

    for (const auto& row : parsingTable) {
        const string& nonTerminal = row.first;
        print << setw(15) << left << nonTerminal;

        for (const string& terminal : terminalList) {
            if (row.second.find(terminal) != row.second.end()) {
                const vector<string>& production = row.second.at(terminal);
                print << setw(15) << left << (nonTerminal + " -> " + production[0]);
            } else {
                print << setw(15) << left << "-";
            }
        }
        print << "\n";
    }
}



//------------------------------- Assignment 3 Logics -------------------------------------------//


// Split a space‑separated line into tokens like ["x","=","5","+",";"]
vector<string> splitTokens(const string &line) {
    vector<string> tokens;
    string current;
    size_t pos = 0;
    while (pos < line.size()) {
        if (isspace(line[pos])) {
            pos++;
            continue;
        }
        // Multi-character terminals (int, if)
        if (pos + 2 < line.size() && line.substr(pos, 3) == "int") {
            tokens.push_back("int");
            pos += 3;
        } else if (pos + 1 < line.size() && line.substr(pos, 2) == "if") {
            tokens.push_back("if");
            pos += 2;
        }
        // Variable (x)
        else if (isalpha(line[pos])) {
            current = line[pos++];
            while (pos < line.size() && isalpha(line[pos])) {
                current += line[pos++];
            }
            tokens.push_back(current);
            current.clear();
        }
        // Number (0, 1, 5)
        else if (isdigit(line[pos])) {
            current = line[pos++];
            tokens.push_back(current);
            current.clear();
        }
        // Single-character terminals (;, =, +, -, >, (, ), {, })
        else {
            tokens.push_back(string(1, line[pos]));
            pos++;
        }
    }
    return tokens;
}

// Print stack contents (bottom -> top)
void printStack(stack<string> st, ostream &out) {
    vector<string> elems;
    while (!st.empty()) {
        elems.push_back(st.top());
        st.pop();
    }
    for (auto it = elems.rbegin(); it != elems.rend(); ++it)
        out << *it << " ";
}

// Breaking a production like "TE'" or "(E)" into ["T","E'"], or ["(","E",")"]
vector<string> tokenizeProduction(const string &prod) {
    vector<string> tokens;
    size_t pos = 0;
    while (pos < prod.size()) {
        if (isspace(prod[pos])) { pos++; continue; }
        // Non‐terminal: uppercase + digits/apostrophe
        if (isupper(prod[pos])) {
            string nt; nt += prod[pos++];
            while (pos < prod.size() && (isdigit(prod[pos]) || prod[pos]=='\'')) {
                nt += prod[pos++];
            }
            tokens.push_back(nt);
        }
        // Terminal: lowercase sequence (e.g. "id")
        else if (islower(prod[pos])) {
            string t;
            while (pos < prod.size() && islower(prod[pos])) {
                t += prod[pos++];
            }
            tokens.push_back(t);
        }
        // Single‐char symbol: + * ( ) etc
        else {
            tokens.push_back(string(1, prod[pos++]));
        }
    }
    return tokens;
}

// Parse one line of input tokens
void parseLine(const vector<string> &inputTokens, const unordered_set<string> &nonTerms, const string &startSymbol, ostream &outConsole, ostream &outFile)
{
    static int lineNo = 1;
    outConsole << "\n--- Parsing line " << lineNo << " ---\n";
    outFile    << "\n--- Parsing line " << lineNo << " ---\n";
    ++lineNo;

    // Initialize stack with [$, startSymbol]
    stack<string> st;
    st.push("$");
    st.push(startSymbol);

    // Prepare input with end‑marker
    vector<string> input = inputTokens;
    input.push_back("$");
    size_t ip = 0;
    int errors = 0;

    // Main parsing loop
    while (!st.empty()) {
        string top = st.top();
        string curr = input[ip];

        // Print step
        outConsole << "Stack: ";  printStack(st, outConsole);
        outConsole << "\tInput: ";
        for (size_t k = ip; k < input.size(); ++k) outConsole << input[k] << " ";
        outConsole << "\nAction: ";

        outFile    << "Stack: ";  printStack(st, outFile);
        outFile    << "\tInput: ";
        for (size_t k = ip; k < input.size(); ++k) outFile << input[k] << " ";
        outFile    << "\nAction: ";

        // 1) Accept when both are "$"
        if (top == "$" && curr == "$") {
            outConsole << "Accept.\n";
            outFile    << "Accept.\n";
            st.pop();
            break;
        }

        // 2) If top is terminal (or "$")
        if (!nonTerms.count(top)) {
            if (top == curr) {
                outConsole << "Match '" << top << "'.\n";
                outFile << "Match '" << top << "'.\n";
                st.pop();
                ip++;
            } else {
                outConsole << "Error: expected '" << top << "', found '" << curr
                           << "'. Skipping input.\n";
                outFile << "Error: expected '" << top << "', found '" << curr
                           << "'. Skipping input.\n";
                errors++;
                if (ip < input.size() - 1) { // Only skip if not at $
                    ip++;
                } else {
                    break; // Stop parsing
                }
            }
        }
        // 3) Top is non‑terminal -> lookup table
        else {
            auto row = parsingTable.find(top);
            bool found = row != parsingTable.end() && row->second.count(curr);
            if (found) {
                // before: const vector<string> &prod = row->second.at(curr);
                // after:
                const vector<string> &rawProd = row->second.at(curr);
                // rawProd should be size 1, holding the whole RHS string
                vector<string> symbols = (rawProd.size()==1)
                    ? tokenizeProduction(rawProd[0])
                    : rawProd;

                // now push symbols in reverse order
                st.pop();
                if (symbols.size()==1 && symbols[0]==EPSILON) {
                    // epsilon‐production
                    outConsole << top << " -> ε\n";
                    outFile    << top << " -> ε\n";
                } 
                else {
                    outConsole << top << " -> ";
                    outFile    << top << " -> ";
                    for (auto &s : symbols) outConsole << s << " ";
                    for (auto &s : symbols) outFile    << s << " ";
                    outConsole << "\n";
                    outFile    << "\n";
                    for (auto it = symbols.rbegin(); it != symbols.rend(); ++it) {
                        st.push(*it);
                    }
                }

            } else {
                outConsole << "Error: no rule for [" << top << ", " << curr
                           << "]. Pop '" << top << "'.\n";
                outFile    << "Error: no rule for [" << top << ", " << curr
                           << "]. Pop '" << top << "'.\n";
                errors++;
                st.pop();
            }
        }
    }

    // Final status
    if (errors==0) {
        outConsole << "Line parsed successfully.\n";
        outFile    << "Line parsed successfully.\n";
    } else {
        outConsole << "Parsing completed with " << errors << " error(s).\n";
        outFile    << "Parsing completed with " << errors << " error(s).\n";
    }
}

// ------------------------------------------------------------------
// Main
int main() {
    string fileName = "input.txt";
    ofstream fileOutput("result.txt");

    vector<Rule> grammar;
    readGrammar(fileName, grammar);
    cout << "\nOriginal Grammar:\n";
    fileOutput << "Original Grammar:\n";
    printGrammar(grammar, cout);
    printGrammar(grammar, fileOutput);

    string originalStartSymbol = grammar[0].lhs;

    leftFactoring(grammar);
    cout << "\nGrammar after Left Factoring:\n";
    fileOutput << "\nGrammar after Left Factoring:\n";
    printGrammar(grammar, cout);
    printGrammar(grammar, fileOutput);

    removeLeftRecursion(grammar);
    cout << "\nGrammar after Left Recursion Removal:\n";
    fileOutput << "\nGrammar after Left Recursion Removal:\n";
    printGrammar(grammar, cout);
    printGrammar(grammar, fileOutput);

    cout << "\nFirst sets of all terminals:\n";
    fileOutput << "\nFirst sets of all terminals:\n";
    computeFirstSets(grammar);
    printFirstSets(cout);
    printFirstSets(fileOutput);

    cout << "\nFollow sets of all terminals:\n";
    fileOutput << "\nFollow sets of all terminals:\n";
    computeFollowSets(grammar, originalStartSymbol);
    printFollowSets(cout);
    printFollowSets(fileOutput);

    cout << "\nParsing table:\n";
    fileOutput << "\nParsing table:\n";
    parsingTable = constructLL1ParsingTable(grammar, firstSet, followSet);
    printParsingTable(parsingTable, cout);
    printParsingTable(parsingTable, fileOutput);

    fileOutput.close();

    unordered_set<string> nonTerms;
    for (auto &r : grammar) nonTerms.insert(r.lhs);

    ifstream in("input_strings.txt");
    if (!in) {
        cerr << "Cannot open input_strings.txt\n";
        return 1;
    }

    ofstream fout("result.txt", ios::app);
    string line;
    vector<string> tokens;
    bool inIfStatement = false;
    while (getline(in, line)) {
        if (trim(line).empty()) continue;
        auto lineTokens = splitTokens(line);
        if (lineTokens.empty()) continue;

        if (!inIfStatement && lineTokens[0] == "if") {
            inIfStatement = true;
            tokens.insert(tokens.end(), lineTokens.begin(), lineTokens.end());
        } else if (inIfStatement) {
            tokens.insert(tokens.end(), lineTokens.begin(), lineTokens.end());
            // Check if line contains '}', assuming single-token '}'
            if (find(lineTokens.begin(), lineTokens.end(), "}") != lineTokens.end()) {
                inIfStatement = false;
                parseLine(tokens, nonTerms, grammar[0].lhs, cout, fout);
                tokens.clear();
            }
        } else {
            parseLine(lineTokens, nonTerms, grammar[0].lhs, cout, fout);
        }
    }

    if (inIfStatement) {
        cerr << "Error: Unclosed if-statement\n";
        fout << "Error: Unclosed if-statement\n";
    }

    return 0;
}