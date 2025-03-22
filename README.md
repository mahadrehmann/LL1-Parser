# LL(1) Parser Construction

As part of our Compiler Construction Assignment 2, we implemented an LL(1) parser by processing a Context-Free Grammar (CFG). It performs the following tasks:

- Reads a CFG from `input.txt`.
- Removes **left recursion**.
- Applies **left factoring**.
- Computes **FIRST** and **FOLLOW** sets.
- Constructs an **LL(1) parsing table**.
- Outputs all results to both the **terminal** and **output.txt**.

## Files

- **LL1-Parsing.cpp** - Main implementation of the LL(1) parser.
- **input.txt** - Contains the CFG in the format: `NonTerminal -> production1 | production2 | ...`.
- **result.txt** - Stores the processing results (grammar transformations, FIRST & FOLLOW sets, parsing table).


## Example Input (input.txt)
```
E -> E+T | T
T -> T*F | F
F -> (E) | id
```
## Output (result.txt)
```
Original Grammar:
E -> E+T | T
T -> T*F | F
F -> (E) | id
-----------------------------------------------------

Grammar after Left Factoring:
E -> E+T | T
T -> T*F | F
F -> (E) | id
-----------------------------------------------------

Grammar after Left Recursion Removal:
E -> TE'
T -> FT'
F -> (E) | id
E' -> +TE' | ^
T' -> *FT' | ^
-----------------------------------------------------

First sets of all terminals:
First(T') = { ^, * }
First(E) = { id, ( }
First(T) = { (, id }
First(F) = { id, ( }
First(E') = { ^, + }
-----------------------------------------------------

Follow sets of all terminals:
Follow(T') = { +, ), $ }
Follow(E) = { ), $ }
Follow(T) = { ), $, + }
Follow(F) = { +, ), $, * }
Follow(E') = { ), $ }
-----------------------------------------------------

Parsing table:
Non-Terminal   $              (              )              *              +              id             
---------------------------------------------------------------------------------------------------------
T'             T' -> ^        -              T' -> ^        T' -> *FT'     T' -> ^        -              
E              -              E -> TE'       -              -              -              E -> TE'       
T              -              T -> FT'       -              -              -              T -> FT'       
F              -              F -> (E)       -              -              -              F -> id        
E'             E' -> ^        -              E' -> ^        -              E' -> +TE'     -              
```
