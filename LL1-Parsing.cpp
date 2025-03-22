#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <algorithm>

using namespace std;

const string EPSILON = "^";
unordered_map<string, unordered_set<string>> firstSet;
unordered_map<string, unordered_set<string>> followSet;

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
void printGrammar(const vector<Rule> &grammar) {
    cout << "-----------------------------------------------------\n";
    for (const auto &rule : grammar) {
        cout << rule.lhs << " -> ";
        for (size_t i = 0; i < rule.productions.size(); i++) {
            cout << rule.productions[i];
            if(i != rule.productions.size() - 1)
                cout << " | ";
        }
        cout << endl;
    }
    cout << "-----------------------------------------------------\n";
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
        vector<string> recursivePart;       // recursive parts
        vector<string> nonRecursivePart;   // non-recursive parts

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

            string newNonTerminal = grammar[i].lhs + "'";
            vector<string> updatedProduction;

            for (auto &betas : nonRecursivePart) {
                if (betas == EPSILON)
                    updatedProduction.push_back(newNonTerminal);
                else
                    updatedProduction.push_back(betas + newNonTerminal);
            }

            grammar[i].productions = updatedProduction;

            Rule newRule;  //Full new rule
            newRule.lhs = newNonTerminal;

            for (auto &alphas : recursivePart) {
                if (alphas == EPSILON)
                    newRule.productions.push_back(newNonTerminal);
                else
                    newRule.productions.push_back(alphas + newNonTerminal);
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
void printFirstSets() {
    cout << "-----------------------------------------------------\n";
    for (const auto &entry : firstSet) {
        cout << "First(" << entry.first << ") = { ";
        bool first = true;
        for (const string &sym : entry.second) {
            if (!first) cout << ", ";
            cout << sym;
            first = false;
        }
        cout << " }" << endl;
    }
    cout << "-----------------------------------------------------\n";
}

// Computing follow sets for all non terminals
void computeFollowSets(vector<Rule> &grammar) {

    for (const auto &rule : grammar) {
        followSet[rule.lhs] = {};
    }
    if (!grammar.empty())
        followSet[grammar[0].lhs].insert("$");


    bool changed;
    do {
        changed = false;

        for (const auto &rule : grammar) {
            string A = rule.lhs;
          
            for (const string &prod : rule.productions) {

                vector<string> symbols;
                size_t pos = 0;
                while (pos < prod.size()) {
                    string token;

                    if (isupper(prod[pos])) {
                        token.push_back(prod[pos++]);
                        if (pos < prod.size() && prod[pos] == '\'') {
                            token.push_back(prod[pos++]);
                        }
                        symbols.push_back(token);
                    }
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
                        
                        // first nikalne ki logic
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
void printFollowSets() {
    cout << "-----------------------------------------------------\n";
    for (const auto &entry : followSet) {
        cout << "Follow(" << entry.first << ") = { ";
        bool first = true;
        for (const string &sym : entry.second) {
            if (!first) cout << ", ";
            cout << sym;
            first = false;
        }
        cout << " }" << endl;
    }
    cout << "-----------------------------------------------------\n";
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
void printParsingTable(const unordered_map<string, unordered_map<string, vector<string>>>& parsingTable) {
    unordered_set<string> terminals;

    for (const auto& row : parsingTable) {
        for (const auto& entry : row.second) {
            terminals.insert(entry.first);
        }
    }

    vector<string> terminalList(terminals.begin(), terminals.end());
    sort(terminalList.begin(), terminalList.end());

    cout << setw(15) << left << "Non-Terminal";
    for (const string& terminal : terminalList) {
        cout << setw(15) << left << terminal;
    }
    cout << "\n" << string(15 + 15 * terminalList.size(), '-') << "\n";

    for (const auto& row : parsingTable) {
        const string& nonTerminal = row.first;
        cout << setw(15) << left << nonTerminal;

        for (const string& terminal : terminalList) {
            if (row.second.find(terminal) != row.second.end()) {
                const vector<string>& production = row.second.at(terminal);
                cout << setw(15) << left << (nonTerminal + " -> " + production[0]);
            } else {
                cout << setw(15) << left << "-";
            }
        }
        cout << "\n";
    }
}

// Main function.
int main() {
    string fileName = "input.txt";
    vector<Rule> grammar;
    unordered_map<string, unordered_map<string, vector<string>>> parsingTable;

    readGrammar(fileName, grammar);
    cout << "Original Grammar:\n";
    printGrammar(grammar);

    // Step 2: Left Factoring
    leftFactoring(grammar);
    cout << "\nGrammar after Left Factoring:\n";
    printGrammar(grammar);

    // Step 3: Left Recursion Removal
    removeLeftRecursion(grammar);
    cout << "\nGrammar after Left Recursion Removal:\n";
    printGrammar(grammar);

    // Step 4: First Sets
    cout << "\nFirst sets of all terminals:\n";
    computeFirstSets(grammar);
    printFirstSets();

    // Step 5: FOLLOW Sets
    cout << "\nFollow sets of all terminals:\n";
    computeFollowSets(grammar);
    printFollowSets();
    
    // Step 6: LL(1) Parsing Table
    parsingTable = constructLL1ParsingTable(grammar, firstSet, followSet);
    printParsingTable(parsingTable);

    return 0;
}
